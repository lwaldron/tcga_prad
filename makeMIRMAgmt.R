######  Converting Refseq identifiers to HGNC symbols with biomaRt ###########

library(readr)

download.file("http://mirdb.org/miRDB/download/miRDB_v5.0_prediction_result.txt.gz", destfile = "miRDB_v5.0_prediction_result.txt.gz")
con <- gzfile("miRDB_v5.0_prediction_result.txt.gz")
dbdat <- readr::read_delim(con, delim = "\t", col_names = FALSE)
dbdat <- data.frame(dbdat, stringsAsFactors = FALSE)

# mart = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
# listFilters(mart)
# ann <- getBM(attributes = c("refseq_mrna", "hgnc_symbol"), filters = "refseq_mrna", values = unique(dbdat[, 2]), mart=mart)
# summary(unique(dbdat[, 2]) %in% ann[, 1])

library(org.Hs.eg.db)
##keytypes(org.Hs.eg.db)
##columns(org.Hs.eg.db)
map <- mapIds(org.Hs.eg.db, keys=dbdat[, 2], column=c("SYMBOL"), keytype="REFSEQ")
dbdat$hgnc <- map
dbdat <- dbdat[!is.na(dbdat$hgnc), ]

mirs <- unique(dbdat[, 1])

mirgmt <- lapply(mirs, function(mir){
  return( dbdat[dbdat[, 1] %in% mir, "hgnc"] )
})
names(mirgmt) <- mirs

writeGMT <- function (object, fname){
  if (class(object) != "list") 
    stop("object should be of class 'list'")
  if (file.exists(fname)) 
    unlink(fname)
  for (iElement in 1:length(object)) {
    write.table(t(c(rep(names(object)[iElement], 2), object[[iElement]])), 
                sep = "\t", quote = FALSE, file = fname, append = TRUE, 
                col.names = FALSE, row.names = FALSE)
  }
}

writeGMT(mirgmt, fname = "mirDB.gmt")

