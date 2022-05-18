
// before running this script, you need to manually concatenate the demultiplexed fastq files. 
// This script expects one fastq file per genome. guppy can do the barcode demultiplexing during basecalling
// now and creates lots of small fastq (or fastq.gz) files in folders called 'barcodeXX', where
// 'XX' stands for the barcode number i.e. 'barcode06'
// For this script to work and name everything correctly you need to concatenate all those files into one .fastq 
// file, NOT fastq.gz!
// There might be a problem of read duplication. To be sure run this bit of code over the concatenated fastq files
//
// for sample in `ls *.fastq | cut -f1 -d'.'`; do cat $sample.fastq | seqkit rmdup -n -o $sample.clean.fastq; done

def helpMessage() {
    log.info"""
    # Nanopore genome polishing
    A pipeline for polishing genomes assembled from Oxford Nanopore reads using Racon and Medaka.

    ## Examples
    nextflow run nanopore_polishing.nf \
    --nanoporeReads "03-trimmed-fastq/*.fastq.gz"
    

    ## Parameters
    --nanoporeReads <glob>
        Required
        A glob of the fastq.gz files of the adapter and barcode trimmed reads.
        The basename of the file needs to match the basename of the respective genome.

    --outdir <path>
        Default: `assembly`
        The directory to store the results in.

    ## Exit codes
    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}

if ( params.nanoporeReads ) {
    nanoporeReads = Channel
    .fromPath(params.nanoporeReads, checkIfExists: true, type: "file")
    .map {file -> [file.simpleName, file]}
    .tap { readsForAssembly }
} else {
    log.info "No nanopore reads supplied, use '--nanoporeReads' and make sure to include '*.fastq.gz'"
    exit 1
}


process canu_version {

    label "canu"
    container "quay.io/biocontainers/canu:2.2-0"

    tag {sampleID} 

    output:
    file 'versions.txt' into version_canu

    """
    echo canu: >> versions.txt
    canu --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process minimap2_version {

    label "minimap2"
    container "quay.io/biocontainers/minimap2:2.24-1"

    tag {sampleID}

    output:
    file 'versions.txt' into version_minimap2

    """
    echo minimap2: >> versions.txt
    minimap2 --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process racon_version {

    label "racon"
    container "quay.io/biocontainers/medaka:1.6.0-0"

    tag {sampleID}

    output:
    path 'versions.txt' into version_racon

    """
    echo racon: >> versions.txt
    racon --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process medaka_version {

    container "quay.io/biocontainers/medaka:1.6.0-0"

    tag {sampleID}

    output:
    path 'versions.txt' into medaka_version

    """
    echo medaka: >> versions.txt
    medaka --version >> versions.txt
    echo --------------- >> versions.txt
    """
}

process seqkit_version {

    container "quay.io/biocontainers/seqkit:2.2.0-0"

    tag {sampleID}

    output:
    file 'versions.txt' into seqkit_version

    """
    echo seqkit: >> versions.txt
    seqkit version >> versions.txt
    echo --------------- >> versions.txt
    """
}


process version {

    input:
    path "canu" from canu_version
    path "racon" from racon_version
    path "minimap" from minimap2_version
    path "medaka" from medaka_version
    path "seqkit" from seqkit_version

    publishDir "${params.outdir}/", mode: 'copy', pattern: 'versions.txt'

    output:
    path "versions.txt"

    script:
    """
    cat canu racon medaka minimap seqkit > versions.txt
    """
}

// genome assembly
process canu {

    container "quay.io/biocontainers/canu:2.2-0"

    tag {sampleID}
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta'
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.fasta.gz'
    publishDir "${params.outdir}/04-canu-assembly", mode: 'copy', pattern: '*.report'

    memory '30 GB'

    input:
    set sampleID, 'input.fastq.gz' from readsForAssembly

    output:
    set sampleID, "${sampleID}.contigs.fasta", 'input.fastq.gz' into minimap2
    set sampleID, "${sampleID}.correctedReads.nanopore.fasta.gz" into correctedReads
    file "${sampleID}.canu.nanopore.report"

    """
    canu \
    -p ${sampleID} \
    -d ${sampleID} \
    genomeSize=45m \
    minInputCoverage=5 \
    stopOnLowCoverage=5 \
    -fast \
    -nanopore input.fastq.gz

    cp ${sampleID}/*contigs.fasta ${sampleID}.contigs.fasta
    cp ${sampleID}/*correctedReads.fasta.gz ${sampleID}.correctedReads.nanopore.fasta.gz
    cp ${sampleID}/*.report ${sampleID}.canu.nanopore.report
    """

}

process minimap2 {

    tag {sampleID}
    label "minimap2"
    container "quay.io/biocontainers/minimap2:2.24-1"

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from minimap2

    output:
    set sampleID, 'input.fasta', 'input.fastq.gz', 'minimap.racon.paf' into racon

    """
    minimap2 \
        input.fasta \
        input.fastq.gz > minimap.racon.paf
    """
}

// polishing step 1
process racon {

    container "quay.io/biocontainers/racon:1.5.0-0"

    tag {sampleID}
    publishDir "${params.outdir}/05-racon-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz', 'minimap.racon.paf' from racon

    output:
    set sampleID, "${sampleID}.contigs.racon.fasta", 'input.fastq.gz' into medaka

    """
    racon -m 8 -x -6 -g -8 -w 500 -t 14\
    --no-trimming \
    input.fastq.gz \
    minimap.racon.paf \
    input.fasta > ${sampleID}.contigs.racon.fasta
    """
}

// polishing step 2
process medaka {

    container "quay.io/biocontainers/medaka:1.6.0-0"

    tag {sampleID} 
    publishDir "${params.outdir}/06-medaka-polish", mode: 'copy', pattern: '*.fasta'

    input:
    set sampleID, 'input.fasta', 'input.fastq.gz' from medaka

    output:
    set sampleID, "${sampleID}.contigs.racon.medaka.fasta", 'input.fastq.gz' into pilon

    """
    medaka_consensus \
    -d input.fasta \
    -i input.fastq.gz \
    -o ${sampleID}_medaka_output \
    -m r941_min_sup_g507

    seqkit sort -lr ${sampleID}_medaka_output/consensus.fasta > ${sampleID}.contigs.racon.medaka.fasta
    seqkit replace -p '.+' -r '${sampleID}_ctg_{nr}' --nr-width 2 ${sampleID}.fasta > ${sampleID}.contigs.racon.medaka.fasta
    """
}
