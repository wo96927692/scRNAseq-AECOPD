samples = ["H4", "H5", "H6", "H2", "H3", "H10", "F11", "H1", "G6", "G10", "G5", "G3", "G4", "G7", "G8", "C24"]
cluster_path = "results/cluster_Myeloid/"
velocity_path = "results/velocity_Myeloid/"
species = "human"
datapath = "/data4/hanzf/run/scs10x/cellranger_results/"
## identify cell markers
celltype_markers = "Immune:CD14"
sub_celltypes = ["Myeloid"]

rule all:
    input: 
        #expand("results/doublet/{sample}_predicted_doublets.txt", sample = samples), 
        #"results/doublet/old/", # update doublet results
        #cluster_path + "seurat_QC.rds", 
        # breakpoint
        cluster_path + "Assess_confounding_factors", 
        cluster_path + "combined_integration.rds", 
        cluster_path + "PCs.csv", 
        cluster_path + "clustree.pdf", 
        cluster_path + "combined_cluster.rds", 
        #"results/ref_markers", 
        #cluster_path + "LabelTransfer", # if the reference seurat data existing
        cluster_path + "markers/all_markers.csv", 
        #cluster_path + "ref_markers_VSL", 
        #cluster_path + "SingleR", # if being mouse data
        cluster_path + "combined_celltypes.rds", 
        directory(cluster_path + "Augur"),
        #directory(cluster_path + "CellChat"), 
        #directory(cluster_path + "CellChat_COMP"), 
        #directory(cluster_path + "celltalker"), 还未完成代码
        #cluster_path + "DEA", 
        cluster_path + "SCENIC/",
        expand("results/cluster_{celltype}/combined_{celltype}.rds", celltype = sub_celltypes), 
        ## velocity
        expand(datapath + "{sample}/velocyto/{sample}.loom", sample = samples),
        directory(velocity_path + "extracted_data"), 
        velocity_path + "myadata_stochastic.h5ad"

rule doublet:
    params:
        datapath = datapath + "{sample}" + "/outs/filtered_feature_bc_matrix",
        samplename = "{sample}"
    output:
        "results/doublet/{sample}_predicted_doublets.txt"
    shell:
        """gunzip -c {params.datapath}/features.tsv.gz > {params.datapath}/features.tsv
           python code/scrublet_doublet.py {params.datapath} {params.samplename} &> logs/doublet.log
        """

rule recallDoublets:
    input:
      expand("results/doublet/{sample}_predicted_doublets.txt", sample = samples)
    params:
        samplename = ",".join(samples)
    output:
        directory("results/doublet/old/")
    shell:
        """Rscript code/recallDoublets.R {params.samplename} 2> logs/recallDoublets.log
        """

rule QC:
    input:
        doublet = expand("results/doublet/{sample}_predicted_doublets.txt", sample = samples)
    params:
        samplename = ",".join(samples), 
        datapath = datapath, 
        results_path = cluster_path
    output:
        cluster_path + "seurat_QC.rds"
    shell:
        """# first argv is samplenames ordered as first col of annotation ("," sperated)
           # second is the path the data stores in (including the last "/")
           Rscript code/QC.R {params.samplename} {params.datapath} {params.results_path}
        """

rule AssessSeurat:
    input:
        #cluster_path + "seurat_QC.rds"
    params:
        results_path = cluster_path, 
        species = species
    output:
        directory(cluster_path + "Assess_confounding_factors"), 
        cluster_path + "seurat_assessed.rds"
    shell:
        """# first argv is samplenames ordered as first col of annotation ("," sperated)
           # second is the path the data stores in (including the last "/")
           Rscript code/AssessSeurat.R {params.results_path} {params.species}
        """

rule integration:
    input:
        after = cluster_path + "Assess_confounding_factors", 
        rds_file = cluster_path + "seurat_assessed.rds"
    params:
        results_path = cluster_path
    output:
        cluster_path + "combined_integration.rds", 
        split_factor = cluster_path + "split_factor.csv"
    shell:
        """# first argv is samplenames ordered as first col of annotation ("," sperated)
           # second is the path the data stores in (including the last "/")
           Rscript code/integration.R {params.results_path} {input.rds_file}
        """

rule PCassessing:
    input: 
        rds_file = cluster_path + "combined_integration.rds", 
        split_factor = cluster_path + "split_factor.csv"
    params:
        results_path = cluster_path
    output:
        cluster_path + "combined_PCassessing.rds", 
        cluster_path + "PCs.csv"
    shell:
        """# first argv is samplenames ordered as first col of annotation ("," sperated)
           # second is the path the data stores in (including the last "/") 
           Rscript code/PCassessing.R {params.results_path} {input.rds_file} {input.split_factor}
        """

rule clustree:
    input:
        rds_file = cluster_path + "combined_PCassessing.rds", 
        PCs = cluster_path + "PCs.csv"
    params:
        results_path = cluster_path
    output:
        cluster_path + "clustree.pdf", 
        cluster_path + "resolution.csv"
    shell:
        """# first argv is numbers of PCs used in download analysis ("," sperated)
           Rscript code/clustree.R {input.rds_file} {params.results_path} {input.PCs}
        """

rule cluster:
    input:
        rds_file = cluster_path + "combined_PCassessing.rds", 
        PCs = cluster_path + "PCs.csv", 
        resolution = cluster_path + "resolution.csv"
    params:
        results_path = cluster_path
    output:
        cluster_path + "combined_cluster.rds"
    shell:
        """# the default resolutions used in finding markers is 1; you have chance to change it
           # the first argv is a python dictionary class including the celltypes and it's markers
           # second argv is numbers of PCs used in download analysis ("," sperated)
           Rscript code/cluster.R {params.results_path} {input.rds_file} {input.PCs} {input.resolution}
        """

rule create_markers:
    input:
    output:
        directory("results/ref_markers")
    shell:
        """Rscript code/conduct_ref_markers.R
        """

rule LabelTransfer:
    input:
        rds_file = cluster_path + "combined_cluster.rds", 
        PCs = cluster_path + "PCs.csv"
    params:
        results_path = cluster_path
    output:
        directory(cluster_path + "LabelTransfer")
    shell:
        """Rscript code/LabelTransfer.R {params.results_path} {input.rds_file} {input.PCs}
        """

rule markers:
    input:
        rds_file = cluster_path + "combined_cluster.rds", 
        PCs = cluster_path + "PCs.csv", 
        resolution = cluster_path + "resolution.csv",
#        afterLT = cluster_path + "LabelTransfer"
    params:
        markers = celltype_markers,
        results_path = cluster_path
    output:
        #cluster_path + "combined_markers.rds", # delete if forget
        cluster_path + "markers/all_markers.csv"
    shell:
        """Rscript code/find_markers.R {params.markers} {params.results_path} {input.rds_file} {input.resolution} {input.PCs}
        """

rule marker_visualization:
    input:
        rds_file = cluster_path + "combined_cluster.rds", 
        afterLT = cluster_path + "LabelTransfer", 
        ref_marker_path = "results/ref_markers", 
        allmarkers_file = cluster_path + "markers/all_markers.csv"
    params:
        results_path = cluster_path
    output:
        directory(cluster_path + "ref_markers_VSL")
    shell:
        """Rscript code/marker_visualization.R {params.results_path} {input.rds_file} {input.ref_marker_path} {input.allmarkers_file}
        """

rule SingleR:
    input:
        rds_file = cluster_path + "combined_cluster.rds",
        afterLT = cluster_path + "LabelTransfer"
    params:
        results_path = cluster_path
    output:
        directory(cluster_path + "SingleR")
    shell:
        """Rscript code/SingleR.R {params.results_path} {input.rds_file} 2> logs/SingleR.log
        """

rule celltypes:
    input: 
        rds_file = cluster_path + "combined_cluster.rds", 
        afterLT = cluster_path + "LabelTransfer",
        markersVSL_path = cluster_path + "ref_markers_VSL"
    params: 
        results_path = cluster_path
    output: 
        cluster_path + "combined_celltypes.rds"
    shell: 
        """Rscript code/celltypes.R {input.rds_file} {params.results_path}
        """

rule Augur: 
    input:
        cluster_path + "combined_celltypes.rds"
    params:
        cluster_path
    output:
        directory(cluster_path + "Augur/")
    shell:
        """Rscript code/Augur.R {params} {input}
	"""

rule CellChat:
    input: 
        rds_file = cluster_path + "combined_celltypes.rds"
    params: 
        results_path = cluster_path, 
        species = species
    output: 
        directory(cluster_path + "CellChat")
    shell: 
        """Rscript code/CellChat.R {params.results_path} {input.rds_file} {params.species}
        """

rule CellChat_comparison:
    input: 
        cluster_path + "CellChat"
    params: 
        results_path = cluster_path, 
    output: 
        directory(cluster_path + "CellChat_COMP")
    shell: 
        """Rscript code/CellChat_comparison.R {params.results_path}
        """

rule celltalker:
    input: 
        rds_file = cluster_path + "combined_celltypes.rds"
    params: 
        species = species
    output: 
        directory(cluster_path + "celltalker")
    shell: 
        """Rscript code/celltalker.R {input.rds_file} {params.species}
        """

rule DE:
    input: 
        rds_file = cluster_path + "combined_celltypes.rds"
    params: 
        results_path = cluster_path
    output: 
        directory(cluster_path + "DEA")
    shell: 
        """Rscript code/DEA.R {input.rds_file} {params.results_path}
        """

rule SCENIC:
    input:
        rds_file = cluster_path + "combined_celltypes.rds"
    params:
        results_path = cluster_path, 
        refdatapath = "/data4/hanzf/run/scs10x/data/SCENIC/", 
        motiffile = "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", 
        motifannota = "motifs-v9-nr.hgnc-m0.001-o0.0.tbl", 
        epigfile = "encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather",
        epigannota = "encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv"
    output:
        directory(protected(cluster_path + "SCENIC/"))
    shell: 
        """Rscript code/SCENIC_makeloomfile.R {input.rds_file} {params.results_path}
           pyscenic grn --num_workers 8 \
                        --output {output}/adj.tsv \
                        --method grnboost2 \
                        {output}/combined_filtered.loom \
                        {params.refdatapath}hs_hgnc_tfs.txt
           pyscenic ctx {output}/adj.tsv \
                    {params.refdatapath}{params.motiffile} \
                    --annotations_fname {params.refdatapath}{params.motifannota} \
                    --expression_mtx_fname {output}/combined_filtered.loom \
                    --mode "dask_multiprocessing" \
                    --output {output}/reg.csv \
                    --num_workers 12 \
                    --mask_dropouts             
           pyscenic aucell \
                    {output}/combined_filtered.loom \
                    {output}/reg.csv \
                    --output {output}/combined_SCENIC.loom \
                    --num_workers 12
        """

rule velocyto:
    input:
        datapath + "{sample}"
    params:
        repeat_masker = "/ssd/public/ref/Hum38_rmsk.gtf", 
        gene_annotation = "/ssd/public/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:
        datapath + "{sample}/velocyto/{sample}.loom"
    shell:
        """velocyto run10x -m {params.repeat_masker} \
           {input} {params.gene_annotation} \
           -@ 4
        """

rule Extract_data:
    input:
        cluster_path + "combined_celltypes.rds"
    params:
        velocity_path
    output:
        directory(velocity_path + "extracted_data"),
        cell_embeddings = velocity_path + "extracted_data/cell_embeddings.csv", 
        cellID_obs = velocity_path + "extracted_data/cellID_obs.csv",
        celltypes = velocity_path + "extracted_data/celltypes.csv", 
        clusters = velocity_path + "extracted_data/clusters.csv"
    shell:
        """Rscript code/Extract_metadata.R {params} {input}
        """

rule scVelo:
    input:
        loom_file = expand(datapath + "{sample}/velocyto/{sample}.loom", sample = samples),
        cell_embeddings = velocity_path + "extracted_data/cell_embeddings.csv", 
        cellID_obs = velocity_path + "extracted_data/cellID_obs.csv",
        celltypes = velocity_path + "extracted_data/celltypes.csv", 
        clusters = velocity_path + "extracted_data/clusters.csv"
    params:
        samples = samples, 
        datapath = datapath, 
        velocity_path = velocity_path
    output:
        velocity_path + "myadata_stochastic.h5ad",
	velocity_path + "myadata_dynamical.h5ad"
    script:
        "scVelo.py"


rule subClusters:
    input: 
        rds_file = ancient(cluster_path + "combined_celltypes.rds")
    params:
        sub_celltype = "{celltype}"
    output:
        "results/cluster_{celltype}/combined_{celltype}.rds"
    shell:
        """Rscript code/subClusters.R {input.rds_file} {params.sub_celltype} {output}
        """
