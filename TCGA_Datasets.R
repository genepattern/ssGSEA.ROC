suppressMessages(suppressWarnings(install.packages("BiocManager", repos = "https://cloud.r-project.org/", 
 quiet = TRUE)))

suppressMessages(suppressWarnings(BiocManager::install("cBioPortalData", ask = FALSE, 
 quiet = TRUE)))
suppressMessages(suppressWarnings(BiocManager::install("AnVIL", ask = FALSE, quiet = TRUE)))
suppressMessages(suppressWarnings(library("cBioPortalData")))
suppressMessages(suppressWarnings(library("AnVIL")))
# suppressMessages(suppressWarnings(install.packages('rapiclient', repos =
# 'https://cloud.r-project.org/', quiet = TRUE))) library('rapiclient')

args = commandArgs(trailingOnly = TRUE)

tcgasamples = as.character(args[2])
symbol_query = as.character(args[4])
threshold_pos = as.numeric(args[6])
threshold_neg = (-1) * abs(as.numeric(args[8]))
msigdbversion = as.character(args[10])

set.seed(147)

cbiosamples <- paste0(tolower(tcgasamples), "_tcga")

# Getting TCGA Dataset
dataset = paste0("http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/", 
 tcgasamples, "/20160128/gdac.broadinstitute.org_", tcgasamples, ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz")
output <- paste0(c("TCGA", tcgasamples, symbol_query, "HIGH_stdev_greater_than", 
 threshold_pos, "vs", "LOW_stdev_less_than_neg", abs(threshold_neg)), collapse = "_")

if (as.numeric(msigdbversion) >= 7.2) {
 symbolchip <- read.table(url(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Human_Gene_Symbol_with_Remapping_MSigDB.v", 
  msigdbversion, ".chip")), header = TRUE, stringsAsFactors = FALSE, sep = "\t", 
  quote = "", fill = TRUE)
} else if (as.numeric(msigdbversion) == 7.1) {
 symbolchip <- read.table(url(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Human_Symbol_with_Remapping_MSigDB.v", 
  msigdbversion, ".chip")), header = TRUE, stringsAsFactors = FALSE, sep = "\t", 
  quote = "", fill = TRUE)
} else if (as.numeric(msigdbversion) < 7.1) {
 message(paste0("Error: MSigDB Version ", msigdbversion, " is not supported. Please try a newer version."))
 stop()
}
symbolchip <- symbolchip[, -c(3)]

symbol_mapped <- symbolchip[symbolchip$Probe.Set.ID == symbol_query, ][, c(2)]
if (length(symbol_mapped) > 1) {
 message("Error: More than one possible mapping was detected for selected gene in the TCGA Dataset")
 print(symbol_mapped)
 stop()
}

temp <- tempfile()
suppressMessages(suppressWarnings(download.file(dataset, temp)))
fname_zipped = basename(dataset)
fnames = as.character(untar(temp, list = TRUE))
untar(temp)
unlink(fnames[basename(fnames) == "MANIFEST.txt"])
fnames = fnames[basename(fnames) != "MANIFEST.txt"]
lst = vector("list", length(fnames))
for (i in seq_along(fnames)) {
 lst[[i]] = read.table(fnames[i], stringsAsFactors = FALSE, sep = "\t")
}
unlink(fnames)
unlink(temp)

data <- lst[[1]]
names <- do.call("rbind", strsplit(data[, c(1)], "|", fixed = TRUE))
data <- cbind(as.data.frame(names[, c(2)], stringsAsFactors = FALSE), data[, data[c(2), 
 ] == "scaled_estimate"])
names(data) = data[c(1), ]
data <- data[-c(1, 2), ]


chip <- read.table(url(paste0("https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Human_NCBI_Entrez_Gene_ID_MSigDB.v", 
 msigdbversion, ".chip")), header = TRUE, stringsAsFactors = FALSE, sep = "\t", 
 quote = "", fill = TRUE)
chip <- chip[, -c(3)]
mappeddata <- merge(x = chip, y = data, by.x = 1, by.y = 1, all = FALSE, no.dups = FALSE)
# ezid_lookup <- mappeddata[mappeddata$Gene.Symbol==symbol_mapped,][,c(1)]
# if(length(ezid_lookup) > 1) { message('More than one possible EntrezGeneID
# mapping was detected for selected gene in the TCGA Dataset') stop() }
mappeddata <- mappeddata[, -c(1)]

mappeddata[, c(2:ncol(mappeddata))] <- sapply(mappeddata[, c(2:ncol(mappeddata))], 
 as.numeric)
mappeddata[is.na(mappeddata)] <- 0

mappeddata <- mappeddata %>% group_by(.data$Gene.Symbol) %>% summarise_all(sum) %>% 
 data.frame()


# Retrieve Sample Expression Thresholding Information from cBioPortal

cbiodata <- suppressMessages(suppressWarnings(cBioDataPack(cbiosamples, ask = FALSE)))
assays <- assays(cbiodata)
cbioassay <- assays$RNA_Seq_v2_mRNA_median_all_sample_Zscores
cbioassay <- as.data.frame(cbind(rownames(cbioassay), cbioassay), stringsAsFactors = FALSE)
mappedcbioassay <- merge(x = symbolchip, y = cbioassay, by.x = 1, by.y = 1, all = FALSE, 
 no.dups = FALSE)
mappedcbioassay <- mappedcbioassay[, -c(1)]
mappedcbioassay[, c(2:ncol(mappedcbioassay))] <- sapply(mappedcbioassay[, c(2:ncol(mappedcbioassay))], 
 as.numeric)
mappedcbioassay[is.na(mappedcbioassay)] <- 0
mappedcbioassay <- mappedcbioassay %>% group_by(.data$Gene.Symbol) %>% summarise_all(max) %>% 
 data.frame()
allnames <- names(mappeddata)

sample_pos <- names(mappedcbioassay)[mappedcbioassay[mappedcbioassay$Gene.Symbol == 
 symbol_mapped, ] >= as.numeric(threshold_pos)]
sample_pos <- sample_pos[sample_pos != "Gene.Symbol"]
sample_neg <- names(mappedcbioassay)[mappedcbioassay[mappedcbioassay$Gene.Symbol == 
 symbol_mapped, ] <= as.numeric(threshold_neg)]
sample_neg <- sample_neg[sample_neg != "Gene.Symbol"]


matches_pos <- unique(grep(paste(sample_pos, collapse = "|"), allnames, value = TRUE))
matches_neg <- unique(grep(paste(sample_neg, collapse = "|"), allnames, value = TRUE))
fullnames <- c(matches_pos, matches_neg)
restricted <- mappeddata[, fullnames]

clslabels <- c(replicate(length(matches_pos), "HIGH"), replicate(length(matches_neg), 
 "LOW"))

write.table(paste(c(dim(restricted)[2], 2, "1"), collapse = " "), paste0(output, 
 ".cls"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(paste("#", paste0(symbol_query, "_HIGH"), paste0(symbol_query, "_LOW"), 
 collapse = " "), paste0(output, ".cls"), sep = "\t", row.names = FALSE, col.names = FALSE, 
 quote = FALSE, append = TRUE)
write.table(paste(clslabels, collapse = " "), paste0(output, ".cls"), sep = "\t", 
 row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

matrix <- as.data.frame(cbind(NAME = mappeddata[, c(1)], Description = mappeddata[, 
 c(1)], restricted), stringsAsFactors = FALSE)

write.table("#1.2", paste0(output, ".TPM.gct"), row.names = FALSE, col.names = FALSE, 
 quote = FALSE)
write.table(t(as.data.frame(dim(restricted))), paste0(output, ".TPM.gct"), sep = "\t", 
 row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
suppressWarnings(write.table(matrix, paste0(output, ".TPM.gct"), sep = "\t", row.names = FALSE, 
 col.names = TRUE, quote = FALSE, append = TRUE))
