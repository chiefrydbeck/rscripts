##load libraries
##to work wiht heatmap
library(gplots)
library(ape)

##Read blast data for all 67 fishes agians PXR
completeFishSet <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/05_17Jun_2016/completeSet/completeSet.txt", header=TRUE, sep="\t", row.names = 1)#, row.names = 1

##read in name key between CEES specific latin ´, american and norwegain names etc
fishNamesCompleteSet <- read.table(file="/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/05_17Jun_2016/completeSet/fishNames_3_completeSet.txt", header=TRUE, sep="\t")#, row.names = 1

##create color vector for heatmap
rgb.paletteYellowGreen <- colorRampPalette(c("grey25", "green"), space = "rgb") 
rgb.paletteRedYellow <- colorRampPalette(c("red", "grey25"), space = "rgb") 
colGradYellowGreen <- rgb.paletteYellowGreen(34)
colGradRedYellow <- rgb.paletteRedYellow(1)
colVec2 <- c(colGradRedYellow,colGradYellowGreen)#colGradRedYellow,

####
##read in tree structure
fishTree <- read.nexus("/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/05_17Jun_2016/completeSet/tree/76g_nucl_conc_unconst.combined.tre")

##https://www.biostars.org/p/134182/
##turn the phylo tree to a dendrogram object
hc <- as.hclust(fishTree) #Compulsory step as as.dendrogram doesn't have a method for phylo objects.
dend <- as.dendrogram(hc)

mergeTableCompleteFishSet <- merge(fishNamesCompleteSet,completeFishSet,by=1)
row.names(mergeTableCompleteFishSet) <- mergeTableCompleteFishSet[,4]

##insert row names before first column of df so that it can be merged with "names df"
completeFishSet_2 <- cbind(completeFishSet[,0:0,drop=F], row.names(completeFishSet), completeFishSet[,(0):length(completeFishSet),drop=F])

##insert other than CEES name by merge names
mergeTableCompleteFishSet_2 <- merge(fishNamesCompleteSet,completeFishSet_2,by=1)

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

##Use latin names (instead of CEES) for row names
row.names(combined_ordered_matrix) <- combined_ordered_matrix[,4]

##key.xtickfun
##function computing tick location and labels for the xaxis of the color key. Returns a named list containing parameters that can be passed to axis. See examples.

#pdf(file = "/Users/halfdanr/NSC/06_Service_Downstream/03_AtlanticCodChemicalDefensome/completeSet/completeSet_withDendrogram_test.pdf",height=5, width=2,pointsize=0.8)
heatmap.2(as.matrix(combined_ordered_matrix[,7:length(combined_ordered_matrix)]), dendrogram = "row", col=colVec2, main="Blast hits per \n position in PXR", trace=c("none"), Rowv=dend, density.info=c("none"), Colv = FALSE, cexCol=0.2,cexRow=0.45,  labCol = "", xlab = "AA position 1-430 in PXR", ,keysize=2) #
#dev.off()
