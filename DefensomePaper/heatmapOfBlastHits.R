#pxr_scf_parsTable <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/66_fishes/70_PXR_scf_parseTable_mod2.txt", header=TRUE, sep="\t")#, row.names = 1
#pxr_scf_parsTable_2 <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/66_fishes/70_PXR_scf_parseTable_mod2.txt", header=TRUE, sep="\t", row.names = 1)#
##load libraries
##to work wiht heatmap
library(gplots)
#to work with dendrograms
library(ape)
#library(phylobase)

##Read blast data for all 67 fishes agians PXR
completeFishSet <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/completeSet.txt", header=TRUE, sep="\t", row.names = 1)#, row.names = 1

#fishNames <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/66_fishes/fishNames_3.txt", header=TRUE, sep="\t")#, row.names = 1
##read in name key between CEES specific latin ´, american and norwegain names etc
fishNamesCompleteSet <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/fishNames_3_completeSet.txt", header=TRUE, sep="\t")#, row.names = 1

#head(fishNames)
#head(pxr_scf_parsTable)


#row.names(pxr_scf_parsTable)

#fishNames[,1]

#dim(fishNames)

#mergeTable <- merge(fishNames,pxr_scf_parsTable,by=1)

#head (mergeTable)
##create color vector for heatmap
rgb.paletteYellowGreen <- colorRampPalette(c("grey25", "green"), space = "rgb") 
rgb.paletteRedYellow <- colorRampPalette(c("red", "grey25"), space = "rgb") 
colGradYellowGreen <- rgb.paletteYellowGreen(34)
colGradRedYellow <- rgb.paletteRedYellow(1)
colVec2 <- c(colGradRedYellow,colGradYellowGreen)#colGradRedYellow,

# Row clustering
# Complete
#pearson corr as dist
#hr.complete.pearson <- hclust(as.dist(1-cor(t(as.matrix(exp)), method="pearson")), method="complete") # Clusters rows by 		Pearson correlation 	as distance method. Pearson better when relationship is linear
#spearman corr as dist
#hr.complete.spearman <- hclust(as.dist(1-cor(t(as.matrix(exp)), method="spearman")), method="complete")
#euclidian dist
#hr.complete.euclidean <- hclust(dist(t(as.matrix(mergeTable)), method="euclidean"), method="complete") 

#rownames <- c(fishNames[,3])
#row.names(mergeTable) <- mergeTable[,4]



#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/66_fishes/hm.pdf",height=5, width=2,pointsize=0.8)
#heatmap.2(as.matrix(mergeTable[,7:436]), dendrogram = "row", col=colVec2, main="Blast hits \n per aa position", trace=c("none"), density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45) #
#dev.off()

#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/66_fishes/hm_2.pdf",height=5, width=2,pointsize=0.8)
#heatmap.2(as.matrix(pxr_scf_parsTable_2), dendrogram = "row", col=colVec2, main="Blast hits \n per aa position", trace=c("none"), density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45) #
#dev.off()

#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/completeSet.pdf",height=5, width=2,pointsize=0.8)
#heatmap.2(as.matrix(completeFishSet), dendrogram = "row", col=colVec2, main="Blast hits \n per aa position", trace=c("none"), density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45) #
#dev.off()

####
##read in tree structure
fishTree <- read.nexus("/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/tree/76g_nucl_conc_unconst.combined.tre")

#plot(fishTree)

##https://www.biostars.org/p/134182/
##turn the phylo tree to a dendrogram object
hc <- as.hclust(fishTree) #Compulsory step as as.dendrogram doesn't have a method for phylo objects.
dend <- as.dendrogram(hc)
#plot(dend, horiz=TRUE) #check dendrogram face

#check Dendrogram
#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/dendrogram.pdf",height=5, width=2,pointsize=0.8)
#plot(dend, horiz=TRUE) #check dendrogram face
#dev.off()

mergeTableCompleteFishSet <- merge(fishNamesCompleteSet,completeFishSet,by=1)
row.names(mergeTableCompleteFishSet) <- mergeTableCompleteFishSet[,4]

#http://www.polarmicrobes.org/?p=562
#rep_tree_r <- root(rep_tree,
#                   resolve.root = T,
#                   interactive = T
#)


##insert row names before first column of df so that it can be merged with "names df"
completeFishSet_2 <- cbind(completeFishSet[,0:0,drop=F], row.names(completeFishSet), completeFishSet[,(0):length(completeFishSet),drop=F])
#head(completeFishSet_2)
##insert other than CEES name by merge names
mergeTableCompleteFishSet_2 <- merge(fishNamesCompleteSet,completeFishSet_2,by=1)
#head(mergeTableCompleteFishSet_2)

##Need to set row names from numbered to CEES Names (Did row names dissappear at merge?)
row.names(mergeTableCompleteFishSet_2) <- mergeTableCompleteFishSet_2[,1]

##http://www.polarmicrobes.org/?p=562
##force row order so that it matches the order of leafs in rep_tree_d
clade_order <- order.dendrogram(dend)
clade_name <- labels(dend)
clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]
new_order <- match(clade_position$clade_name, row.names(mergeTableCompleteFishSet_2))
combined_ordered_matrix <- mergeTableCompleteFishSet_2[new_order,]
#head(combined_ordered_matrix)

##Use latin names (instead of CEES) for row names
row.names(combined_ordered_matrix) <- combined_ordered_matrix[,4]


#head (combined_ordered_matrix)


#row.names(combined_ordered_matrix)
##key.xtickfun
##function computing tick location and labels for the xaxis of the color key. Returns a named list containing parameters that can be passed to axis. See examples.

pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/completeSet_withDendrogram_test.pdf",height=5, width=2,pointsize=0.8)
heatmap.2(as.matrix(combined_ordered_matrix[,7:length(combined_ordered_matrix)]), dendrogram = "row", col=colVec2, main="Blast hits per \n position in PXR", trace=c("none"), Rowv=dend, density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45,  labCol = "", xlab = "AA position 1-430 in PXR", ,keysize=2) #
dev.off()

#"Fish species  \n blast hits per  \n amino acid position \n in PXR"


#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/completeSet_withDendrogram.pdf",height=5, width=2,pointsize=0.8)
#heatmap.2(as.matrix(mergeTableCompleteFishSet[,7:436]), dendrogram = "row", col=colVec2, main="Blast hits \n per aa position", trace=c("none"), Rowv=dend, density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45) #
#dev.off()


#dim(mergeTableCompleteFishSet)
#dim(merge(fishNamesCompleteSet[,1],completeFishSet[,1], by=1))
