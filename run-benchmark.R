##########
## Info ##
##########

# The code was run on a computer with the following setup:
# 64 GB RAM
# 32 cores Intel Core Processor (Broadwell) 2.6 GHz
# Ubuntu 16.04.3 64 bit
# Python 3.6.4
# R 3.4.3
#
# All methods are up to date on 20/03/2018.
# The output of sessionInfo() in R can be found at the bottom of this document

################################
## Generate Artifical Dataset ##
################################

# generate artifical dataset from hematopoisis dataset by Velten et al.
library(splatter)

cmat <- data.matrix(read.csv("data/GSE75478_transcriptomics_raw_filtered_I1.csv.gz", row.names=1))
params <- splatEstimate(cmat)
saveRDS(params, "data/splatter-params.rds")

###################
## Run Benchmark ##
###################

source("code/benchmark-methods.R")

# cell scaling
seeds <- c(464, 765, 694, 902, 111, 472, 284, 492, 534, 889)
nCellsVec <- c(50, 200, 500, 1000, 2000, 5000, 10000, 20000, 30000)
runBatch("slingshot", seeds, 5000, nCellsVec)
runBatch("TSCAN", seeds, 5000, nCellsVec)
runBatch("DTP", seeds, 5000, nCellsVec)
runBatch("monocle", seeds, 5000, nCellsVec)
runBatch("scanpyAGA", seeds, 5000, nCellsVec)
runBatch("scanpyDPT", seeds, 5000, nCellsVec)
runBatch("wishbone", seeds, 5000, nCellsVec)
runBatch("GPfates", seeds, 5000, nCellsVec)


# gene scaling
seeds <- c(464, 765, 694, 902, 111, 472, 284, 492, 534, 889)
nGenesVec <- c(50, 100, 500, 1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000)
runBatch("slingshot", seeds, nGenesVec, 500)
runBatch("TSCAN", seeds, nGenesVec, 500)
runBatch("DTP", seeds, nGenesVec, 500)
runBatch("monocle", seeds, nGenesVec, 500)
runBatch("scanpyAGA", seeds, nGenesVec, 500)
runBatch("scanpyDPT", seeds, nGenesVec, 500)
runBatch("wishbone", seeds, nGenesVec, 500)
runBatch("GPfates", seeds, nGenesVec, 500)




# sessionInfo()
# R version 3.4.3 (2017-11-30)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
#
# Matrix products: default
# BLAS: /usr/lib/libblas/libblas.so.3.6.0
# LAPACK: /usr/lib/lapack/liblapack.so.3.6.0
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#
# attached base packages:
#  [1] splines   stats4    parallel  stats     graphics  grDevices utils
#  [8] datasets  methods   base
#
# other attached packages:
#  [1] destiny_2.6.2              TSCAN_1.16.0
#  [3] monocle_2.6.3              DDRTree_0.1.5
#  [5] irlba_2.3.2                VGAM_1.0-5
#  [7] Matrix_1.2-12              slingshot_0.1.2-3
#  [9] princurve_1.1-12           tictoc_1.0
# [11] reticulate_1.6             splatter_1.2.1
# [13] scater_1.6.3               SingleCellExperiment_1.0.0
# [15] SummarizedExperiment_1.8.1 DelayedArray_0.4.1
# [17] matrixStats_0.53.1         GenomicRanges_1.30.3
# [19] GenomeInfoDb_1.14.0        IRanges_2.12.0
# [21] S4Vectors_0.16.0           ggplot2_2.2.1
# [23] Biobase_2.38.0             BiocGenerics_0.24.0
#
# loaded via a namespace (and not attached):
#   [1] backports_1.1.2        Hmisc_4.1-1            RcppEigen_0.3.3.4.0
#   [4] plyr_1.8.4             igraph_1.2.1           lazyeval_0.2.1
#   [7] sp_1.2-7               shinydashboard_0.7.0   BiocParallel_1.12.0
#  [10] densityClust_0.3       fastICA_1.2-1          digest_0.6.15
#  [13] htmltools_0.3.6        viridis_0.5.0          gdata_2.18.0
#  [16] magrittr_1.5           checkmate_1.8.5        memoise_1.1.0
#  [19] cluster_2.0.6          limma_3.34.9           xts_0.10-2
#  [22] prettyunits_1.0.2      colorspace_1.3-2       blob_1.1.0
#  [25] ggrepel_0.7.0          dplyr_0.7.4            RCurl_1.95-4.10
#  [28] jsonlite_1.5           tximport_1.6.0         lme4_1.1-15
#  [31] bindr_0.1.1            zoo_1.8-1              survival_2.41-3
#  [34] ape_5.0                glue_1.2.0             gtable_0.2.0
#  [37] zlibbioc_1.24.0        XVector_0.18.0         MatrixModels_0.4-1
#  [40] car_2.1-6              DEoptimR_1.0-8         SparseM_1.77
#  [43] VIM_4.7.0              scales_0.5.0           pheatmap_1.0.8
#  [46] DBI_0.8                edgeR_3.20.9           Rcpp_0.12.16
#  [49] laeken_0.4.6           viridisLite_0.3.0      xtable_1.8-2
#  [52] progress_1.1.2         htmlTable_1.11.2       proxy_0.4-21
#  [55] foreign_0.8-69         bit_1.1-12             mclust_5.4
#  [58] Formula_1.2-2          vcd_1.4-4              htmlwidgets_1.0
#  [61] httr_1.3.1             FNN_1.1                gplots_3.0.1
#  [64] RColorBrewer_1.1-2     acepack_1.4.1          pkgconfig_2.0.1
#  [67] XML_3.98-1.10          nnet_7.3-12            locfit_1.5-9.1
#  [70] rlang_0.2.0            reshape2_1.4.3         AnnotationDbi_1.40.0
#  [73] munsell_0.4.3          tools_3.4.3            RSQLite_2.0
#  [76] stringr_1.3.0          knitr_1.20             bit64_0.9-7
#  [79] robustbase_0.92-8      caTools_1.17.1         RANN_2.5.1
#  [82] bindrcpp_0.2           nlme_3.1-131.1         quantreg_5.35
#  [85] mime_0.5               slam_0.1-42            biomaRt_2.34.2
#  [88] compiler_3.4.3         pbkrtest_0.4-7         rstudioapi_0.7
#  [91] curl_3.1               beeswarm_0.2.3         e1071_1.6-8
#  [94] smoother_1.1           tibble_1.4.2           stringi_1.1.7
#  [97] lattice_0.20-35        nloptr_1.0.4           HSMMSingleCell_0.112.0
# [100] pillar_1.2.1           lmtest_0.9-35          combinat_0.0-8
# [103] data.table_1.10.4-3    bitops_1.0-6           httpuv_1.3.6.2
# [106] R6_2.2.2               latticeExtra_0.6-28    KernSmooth_2.23-15
# [109] gridExtra_2.3          vipor_0.4.5            boot_1.3-20
# [112] MASS_7.3-49            gtools_3.5.0           assertthat_0.2.0
# [115] rhdf5_2.22.0           rjson_0.2.15           qlcMatrix_0.9.5
# [118] GenomeInfoDbData_1.0.0 mgcv_1.8-23            grid_3.4.3
# [121] rpart_4.1-13           minqa_1.2.4            class_7.3-14
# [124] Rtsne_0.13             TTR_0.23-3             shiny_1.0.5
# [127] base64enc_0.1-3        ggbeeswarm_0.6.0
