####################################################################
# Preprocessing of the trial data (Zhao et al, Nat. Medicine 2019) #
####################################################################

require("xlsx")


## curate annotations

sraAnnot <- fp.quickRead("SraRunTable.csv",sep=",",header=T) ## reading the SRA table for annotation
suppTable <- xlsx::read.xlsx("41591_2019_349_MOESM1_ESM(4).xlsx",1) ## reading from the supplementary table from the paper
authorCorrespondance <- xlsx::read.xlsx("PUBLIC_PD1_TREATED_SAMPLES.xlsx",1) ## sample-matching ID, obtained through correspondance with the authors.

annotTrial <- sraAnnot[sraAnnot$Assay.Type=="RNA-Seq",c(1,2,28,29,33)]
rownames(annotTrial) <- annotTrial$Sample.Name

rownames(authorCorrespondance) <- authorCorrespondance$DNA.Sampel.Name
addAnnot <- authorCorrespondance[rownames(annotTrial),c(1:4,7,10)]
annotTrial <- merge(annotTrial,addAnnot,by="row.names")

colnames(suppTable) <- suppTable[1,]
suppTable <- suppTable[2:nrow(suppTable),]
suppTable <- suppTable[suppTable$`Patient #` %in% as.character(annotTrial$PATIENT_ID),]
suppTable <- suppTable[,c(1,3:5,9:13,16,18:ncol(suppTable))]
annotTrial <- merge(annotTrial,suppTable,by.x="PATIENT_ID",by.y ="Patient #")


rownames(annotTrial) = annotTrial$Run

## gene expression

load("geneCount_raw_31s_.RData") # loads a data.frame called counts with gene counts (ENSG gene IDs)
trial_counts_ENSG <- counts

correspondance <- read.table("20250707_HGNC_symbols.tsv",sep="\t",header = TRUE) ## correspondance table between HUGO gene symbols and ENSEMBL IDs
correspondance <- correspondance[nchar(correspondance$Ensembl.ID.supplied.by.Ensembl.)>3,]
correspondance <- correspondance[nchar(correspondance$Approved.symbol)>1,]
correspondance <- correspondance[!duplicated(correspondance$Ensembl.ID.supplied.by.Ensembl.),]
rownames(correspondance) <- correspondance$Ensembl.ID.supplied.by.Ensembl.

correspondance <- correspondance[intersect(rownames(trial_counts_ENSG),rownames(correspondance)),]
trial_counts_ENSG <- trial_counts_ENSG[rownames(correspondance),]

rownames(trial_counts_ENSG) <- correspondance$Approved.symbol

trial_counts <- trial_counts_ENSG[sort(rownames(trial_counts_ENSG)),]



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



trial_TPM_log2 <- countsTologTPM(trial_counts)
trial_TPM_linear <- (2^trial_TPM_log2)-1


exp_PRJNA482620 <- trial_TPM_log2
annot_PRJNA482620 <- annotTrial

annot_PRJNA482620 <- annot_PRJNA482620[colnames(exp_PRJNA482620),]
rownames(annot_PRJNA482620) <- colnames(exp_PRJNA482620)
