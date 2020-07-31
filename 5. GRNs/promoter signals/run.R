dfs = lapply(1:21, function(i){
    fp = sprintf("D:\\project\\single_cell\\MOB_ATAC\\GRNs\\GRN_data\\promoterSignal_%s.csv", i)
    read.table(fp, header = T, sep = "\t")
})


genes = lapply(dfs, function(df) df$symbol)
