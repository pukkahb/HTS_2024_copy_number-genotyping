#!/usr/bin/env nextflow

nextflow.enable.dsl2

params {
    fastq_dir = 'data/raw/'
    reference = 'data/reference/'
    results_dir = 'results/'
    cnv_dir = 'results/cnv/'
    plots_dir = 'results/plots/'
    threads = 4
}

workflow {
    // Define input channels
    fastq_files = file(params.fastq_dir).glob('*.fastq.gz')
    reference_genome = file(params.reference + 'reference_genome.fa')

    // Quality check
    fastqc_results = file("${params.results_dir}/fastqc")
    multiqc_results = file("${params.results_dir}/multiqc")

    process fastqc_input {
        input:
        path fastq_file from fastq_files

        output:
        path fastqc_results

        script:
        """
        bash scripts/fastqc.sh $fastq_file $fastqc_results $params.threads
        """
    }

    process multiqc {
        input:
        path fastqc_results

        output:
        path multiqc_results

        script:
        """
        bash scripts/multiqc.sh $fastqc_results $multiqc_results ${params.fastq_dir}
        """
    }

    // Trimming
    trimmed_reads = file("${params.results_dir}/trimmed")

    process trim_galore {
        input:
        path fastq_files

        output:
        path trimmed_reads

        script:
        """
        bash scripts/trim_galore.sh $fastq_files $trimmed_reads
        """
    }

    // BWA Indexing
    indexed_reference = file("${params.results_dir}/indexed")

    process bwa_index {
        input:
        path reference_genome

        output:
        path indexed_reference

        script:
        """
        bash scripts/bwa_index.sh $reference_genome $indexed_reference
        """
    }

    // BWA Alignment
    aligned_reads = file("${params.results_dir}/aligned")

    process bwa_align {
        input:
        path trimmed_reads
        path indexed_reference

        output:
        path aligned_reads

        script:
        """
        bash scripts/bwa_align.sh $trimmed_reads $indexed_reference $aligned_reads
        """
    }

    // Variant Calling
    variant_calls = file("${params.results_dir}/variants")

    process samtools_call {
        input:
        path aligned_reads

        output:
        path variant_calls

        script:
        """
        bash scripts/samtools_call.sh $aligned_reads $variant_calls
        """
    }

    // CNV Analysis
    cnv_results = file(params.cnv_dir)

    process cnv_analysis {
        input:
        path variant_calls
        path reference_genome

        output:
        path cnv_results

        script:
        """
        bash scripts/cnv_analysis.sh $variant_calls $reference_genome $cnv_results
        """
    }

    // CNV Plotting
    process plot_cnv {
        input:
        path cnv_results

        output:
        path params.plots_dir

        script:
        """
        Rscript scripts/plot_cnv.R -i $cnv_results -o ${params.plots_dir}
        """
    }
}

