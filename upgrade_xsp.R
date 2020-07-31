# Sometimes you create xsp in SnapATAC-0.9 (old version)
# and you want to use it in SnapATAC-1.0 (new version up
# to 2019.01). You can use the following trick to convert
# old xsp to new object. This trick only reserve the bmat
# of xsp and You have to restart the SnapATAC analysis.

x.sp = readRDS("xsp_0.9.rds")

chroms.v = as.character(x.sp@feature@seqnames@values)
chroms.t = x.sp@feature@seqnames@lengths
chroms = rep(chroms.v, chroms.t)
start = x.sp@feature@ranges@start
bins = GRanges(chroms, IRanges(start, start + 5000 - 1))

xsp = createSnapFromBmat(x.sp@bmat, x.sp@barcode, bins)

xsp@metaData = x.sp@metaData
xsp@file = x.sp@file
xsp@sample = x.sp@sample
xsp@feature = x.sp@feature

rm(list = setdiff(ls(), "xsp"))
saveRDS(xsp, file = "xsp_1.0.rds")