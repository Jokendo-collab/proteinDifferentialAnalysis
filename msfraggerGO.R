
setwd("/home/javan/Desktop/dlovuData/msfragger_results")
dir()
# analysis url: https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/tfxp6w02segq7c6/msfragger_cluser_profiler.nb.html
# http://www.jingege.wang/a/study/r/341.html
# load Libraries (install if necessary) 

 packages = c("BiocManager","tidyverse","clusterProfiler","org.Hs.eg.db")
# 
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# BiocManager::install("clusterProfiler") # may need to uncomment this line  if this does not install because of bioconductor install weirdness 

# BiocManager::install("org.Hs.eg.db") # may need to uncomment this line if this does not install because of bioconductor install weirdness 

#verify they are loaded
search()

# load  msFragger data 
# this is a 100 ng hela run on teh timstofPRO
# Change this to the combined protein file if it's a multiexperiment

data <- read.table("control_protein.tsv",header = T,sep = '\t')

library(org.Hs.eg.db) # this is for human. load a different library from here for a different speccies http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
library(clusterProfiler)
library(enrichplot)

gene <- data$Gene# extract Gene's from MSFragger

# this translates the Gene from MSfragger to something enrichgo can read
gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)


# Make a geneList for some future functions
geneList <- gene.df$ENTREZID
names(geneList) <- as.character(gene.df$SYMBOL)
geneList <- sort(geneList, decreasing = TRUE)


# gene enrichment analysis cnplots are commented out as they look crazy with a large number of proteins
## BP
ego_BP <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   readable = TRUE,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)


ego_BP_df <-as.data.frame(ego_BP) #convert GO terms to dataframe

#Dotplot and barplot only shows the most significantly enriched terms
## remove redundent GO terms first 
ego2 <- simplify(ego_BP)

# Plot either dotplot or barplot
dotplot(ego2, showCategory=20)

barplot(ego2, showCategory=20)

goplot(ego2,foldChange=geneList) #Plot the enriched terms and their significance

# Gene concept network, show the kegg pathways

cnetplot(ego2, foldChange=geneList)

cnetplot(ego2, foldChange=geneList, circular = F, colorEdge = T) #plotting a circular cnetplot, not good with many proteins

#enrich map
# Enrichment map organizes enriched terms into a network with edges connecting overlapping gene sets. 
# In this way, mutually overlapping gene sets are tend to cluster together, 
# making it easy to identify functional module.

emapplot(ego2)


#Cellular component
ego_CC <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   readable = TRUE,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

## remove redundent GO terms
ego3 <- simplify(ego_CC)

ego_CC_df <-as.data.frame(ego_CC)

dotplot(ego3, showCategory=20)

barplot(ego3, showCategory=20)

goplot(ego3,foldChange=geneList) #Plot the enriched terms and their significance
#Plot the cnet to visualize protein cluster
cnetplot(ego3, foldChange=geneList)

cnetplot(ego3, foldChange=geneList, circular = F, colorEdge = TRUE) #plotting a circular cnetplot, not good with many proteins

# Enriched map visualize mutually overlapping proteins
emapplot(ego3)

## MF
ego_MF <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb  = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   readable = TRUE,
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

## remove redundent GO terms
ego4 <- simplify(ego_MF) 

ego_MF_df <-as.data.frame(ego_MF)

dotplot(ego4, showCategory=20)

barplot(ego4, showCategory=20)

goplot(ego4,foldChange=geneList) #Plot the enriched terms and their significance

#Visualize the protein network
cnetplot(ego4, foldChange=geneList, circular = F, colorEdge = TRUE) #plotting a circular cnetplot, not good with many proteins

# Enriched map visualize mutually overlapping proteins
emapplot(ego4)

#Doing all the ontologis at once 
ego_all <- enrichGO(gene = gene.df$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   readable = TRUE,
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

#Combine all the ontology figures together
barplot(ego_all, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

#Combine all the ontology figures together
dotplot(ego_all, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

# Plotting a complex venn diagram using upset plot diagram
upsetplot(ego_BP)

#Biological Process (BP):  
#Cellular Component (CC):
#Molecular Function (MF): 

# # write Go data to csv files
# write_csv(ego_BP_df,"go_biological.process.tbhart.csv")
# write_csv(ego_CC_df,"go_cellular_component.tbhart.csv")
# write_csv(ego_MF_df,"go_molecular_function.tbhart.csv")

# complex venn diagram
upsetplot(ego2)

# Produces heatmap plot
heatplot(ego2)

heatplot(ego2, foldChange=geneList)

#Produces enrichment network plot
emapplot(ego2)

#Ridge plot
gene.df$ENTREZID = bitr(gene.df$SYMBOL, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
geneList = as.numeric(gene.df$ENTREZID)

names(geneList) = as.character(gene.df$ENTREZID)
geneList = sort(geneList, decreasing = T)
#kk2 <- gseKEGG(geneList = geneList,organism = 'hsa', nPerm= 1000, minGSSize = 10, pvalueCutoff = 0.0005,verbose= TRUE)

kk <- gseKEGG(geneList, nPerm=10000)

ridgeplot(kk)


final_productKEGG_Overrep <- heatplot(k1, foldChange = geneList, showCategory 
                                      = 20000)+ ggplot2::coord_flip()



# #=============================================
# library(clusterProfiler)
# 
# library(org.Mm.eg.db)
# data <- UbiLength
# 
# gene <- data$Gene# extract Gene's from MSFragger
# 
# # this translates the Gene from MSfragger to something enrichgo can read
# gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Mm.eg.db)
# 
# #Ridge plot
# gene.df$ENTREZID = bitr(gene.df$SYMBOL, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Mm.eg.db")
# geneList = as.numeric(gene.df$ENTREZID)
# 
# names(geneList) = as.character(gene.df$ENTREZID)
# geneList = sort(geneList, decreasing = T)
# #kk2 <- gseKEGG(geneList = geneList,organism = 'hsa', nPerm= 1000, minGSSize = 10, pvalueCutoff = 0.0005,verbose= TRUE)
# 
# kk <- gseKEGG(geneList, nPerm=10000)
# 
# ridgeplot(kk)



































