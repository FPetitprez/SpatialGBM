############################
# Preparation of TCGA data #
############################

library("SummarizedExperiment")
library("TCGAbiolinks")
library("jsonlite")
library("httr")




# Prepare annotation
annot_TCGA_GBM <- fp.quickRead("gbm_tcga_pub2013/data_clinical_patient.txt",sep="\t") ## read file from CBioPortal
rownames(annot_TCGA_GBM) <- annot_TCGA_GBM$PATIENT_ID

annot_TCGA_GBM$OS.time = annot_TCGA_GBM$OS_MONTHS
annot_TCGA_GBM$OS.event = as.numeric(substr(annot_TCGA_GBM$OS_STATUS,1,1))





# Query platform Illumina HiSeq
query <- GDCquery(
  project = "TCGA-GBM", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

expdat <- GDCprepare(
  query = query,
  save = TRUE,
  save.filename = "TCGA_GBM_exp.rda"
)

TCGA_GBM_counts = assay(expdat,"unstranded") |> as.data.frame()
colnames(TCGA_GBM_counts) = substr(colnames(TCGA_GBM_counts),1,12)


###
# Multiple IDs to convert - use a POST request
###
url = "https://biotools.fr/human/ensembl_symbol_converter/"
ids = rownames(TCGA_GBM_counts)
ids_json <- toJSON(ids)

body <- list(api=1, ids=ids_json)
r <- POST(url, body = body)

output <- fromJSON( content(r, "text"), flatten=TRUE)
geneCorres <- data.frame(ENSG = names(output), Symbol = NA)
for(i in 1:nrow(geneCorres)){
  print(ensg)
  ensg = geneCorres[i,"ENSG"]
  if(!is.null(output[ensg][[1]])) geneCorres[i,"Symbol"] = output[ensg]
}
geneCorres <- geneCorres[!is.na(geneCorres$Symbol),]
rownames(geneCorres) <- geneCorres$ENSG

geneCorres[duplicated(geneCorres$Symbol),"Symbol"] |> sort()

geneCorres <- geneCorres[!duplicated(geneCorres$Symbol),]

TCGA_GBM_counts <- TCGA_GBM_counts[rownames(geneCorres),]
rownames(TCGA_GBM_counts) <- geneCorres$Symbol
TCGA_GBM_counts <- TCGA_GBM_counts[sort(rownames(TCGA_GBM_counts)),]



tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


countsTologTPM =function(counts){
  ## step1: get gene lengths
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
  genes <- rownames(counts)
  G_list <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol", "ensembl_gene_id",
                                                        "description","transcript_length"),values=genes,mart= mart)
  geneLengths = data.frame(geneSymbols = genes, length = NA)
  rownames(geneLengths) = geneLengths$geneSymbols
  for(g in geneLengths$geneSymbols){
    geneLengths[g, "length"] = max(c(subset(G_list,hgnc_symbol == g)$transcript_length),0)
  }

  ## step2: tpm transformation + log scale
  res = counts[rownames(geneLengths),]
  res = tpm3(res, geneLengths$length)
  res = as.data.frame(res)
  res = log2(1+res)
  return(res)
  
}


TCGA_GBM_TPM_log2 <- fp.countsTologTPM(TCGA_GBM_counts)
TCGA_GBM_TPM_linear <- (2^TCGA_GBM_TPM_log2)-1













