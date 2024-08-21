#!/usr/bin/env nextflow

// Initialize parameters with a default value if not provided
params.samples = channel.fromPath( '/data/counts/*.csv' )
params.background_samples = channel.fromPath( '/data/background/*.csv' )
params.fdr = 'fdr'
params.lower = 'lower'
working_dir = params.working_dir ?: "." 


// Create channels from the list of samples
samples_ch = Channel.from(params.samples)
background_samples_ch = Channel.from(params.samples)

// Create channels from all paraters defined in config
fdr_ch = Channel.from(params.fdr)
lower_ch = Channel.from(params.lower)

process Gene_id_to_gene_symbol {
    tag "${sample}"
    
    input:
    val sample

    output:
    path("GeneSymbol_${sample}")

    script:

    """
    echo "adding gene symbol to counts" 
    
    echo "Processing sample: ${sample}"
    
    Rscript ${working_dir}/nextflow_scripts/GeneID_to_geneSymbol.r \\
    --data "${working_dir}/data/counts/${sample}" \\
    --ref_table "${working_dir}/resources/v41_ENSG_ID_GENEsymbol.csv" \\
    --output "GeneSymbol_${sample}"
    
    mkdir -p  ${working_dir}/results/counts 
    cp * ${working_dir}/results/counts  

    """
}

// background process 
process Gene_id_to_gene_symbol_background {
    tag "${sample}"
    
    input:
    val sample

    output:
    //tuple path("background_GeneSymbol_${sample}"), val (sample)
    path("background_GeneSymbol_${sample}")

    script:

    """
    echo "adding gene symbol to background counts" 

    echo "Processing sample: ${sample}"
    
    Rscript ${working_dir}/nextflow_scripts/GeneID_to_geneSymbol.r \\
    --data "${working_dir}/data/background/${sample}" \\
    --ref_table "${working_dir}/resources/v41_ENSG_ID_GENEsymbol.csv" \\
    --output "background_GeneSymbol_${sample}"



    mkdir -p  ${working_dir}/results/background 
    cp * ${working_dir}/results/background 

    """
}

// Process for performing empty drops analysis
process empty_drops_analysis {
    tag "${sample_id}"

    input:
    path gene_symbol_counts 
    path background_gene_symbol_counts
    val sample_id
    val fdr
    val lower

    output:
    path "empty_drops/removed_*"

    
    
    script:
    """
    echo "Running empty drops analysis"
    mkdir -p empty_drops

    Rscript ${working_dir}/nextflow_scripts/empty_drops.r \\
    "empty_drops" \\
    ${gene_symbol_counts} \\
    ${background_gene_symbol_counts} \\
    ${sample_id} \\
    ${fdr} \\
    ${lower}


    mkdir -p  ${working_dir}/results/empty_drops 
    cd empty_drops
    cp * ${working_dir}/results/empty_drops  
    """
}

// Process for performing decontx analysis
process run_decontx {
    tag "${sample_id}"

    input:
    path empty_drops_seurat_obj
    path background_gene_symbol_counts
    val  sample_id

    output:
    path "*_decontx_counts.csv"

    script:
    """
    Rscript ${working_dir}/nextflow_scripts/decontx.r \\
    --seurat_obj "${empty_drops_seurat_obj}" \\
    --background_counts_path "${background_gene_symbol_counts}" \\
    --sample_id "${sample_id}"
    
    mkdir -p ${working_dir}/results/decontx 
    cp * ${working_dir}/results/decontx   
    
    """
}

// Process for performing QC analysis
process gene_QC {
    tag "${sample_id}"

    input:
    path decontx_counts
    val  sample_id

    output:
    stdout

    script:
    """
    Rscript ${working_dir}/nextflow_scripts/QC.R -i "${decontx_counts}" -s "${sample_id}".basename --npc 20 --cluster_res 0.7 --mt 10 --min_counts 1000 --max_counts 100000 --min_features 1000 --max_features 10000000

    mkdir -p ${working_dir}/results/QC/genes
    cp * ${working_dir}/results/QC/genes   
    """
}

// execute workflow
workflow {
    // Execute the processes and collect their outputs
    Gene_id_to_gene_symbol(samples_ch)

    Gene_id_to_gene_symbol_background(background_samples_ch)
    
     // run decontx with four inputs 1. empty drops analysis, 2. background expression 3. sample chanel 
    empty_drops_analysis(Gene_id_to_gene_symbol.out,
                         Gene_id_to_gene_symbol_background.out,
                         samples_ch,
                         fdr_ch, 
                         lower_ch)
    
    // run decontx with three inputs 1. empty drops analysis, 2. background expression 3. sample chanel 
    run_decontx(empty_drops_analysis.out,
                Gene_id_to_gene_symbol_background.out,
                samples_ch)


        // run decontx with three inputs 1. empty drops analysis, 2. background expression 3. sample chanel 
    gene_QC(run_decontx.out,
                samples_ch)
    

    // Assuming `Gene_id_to_gene_symbol_background` outputs directories
    //list_files(background_gene_id).view()
}