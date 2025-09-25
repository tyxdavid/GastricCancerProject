#!/bin/sh

cd ./main_cluster/subcluster/Epithelial/scenic

nohup arboreto_with_multiprocessing.py  seurat_for_monocle.loom \
                    ./cisTarget_databases/hg38/allTFs_hg38.txt \
					-o adjacencies.csv --seed 774 --num_workers 48 \
					 >grn_arboreto.log 2>&1 &

nohup pyscenic ctx \
        adjacencies.csv \
        ./cisTarget_databases/hg38/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        ./cisTarget_databases/hg38/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --annotations_fname ./cisTarget_databases/hg38/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname seurat_for_monocle.loom \
        --mode "custom_multiprocessing" \
        --output regulons.csv \
        --num_workers 48  >ctx.log 2>&1 &

nohup pyscenic aucell \
        seurat_for_monocle.loom \
        regulons.csv \
        -o auc_mtx.csv --seed 774 \
        --num_workers 24  >auc_cell.log 2>&1 &
