##Master's Thesis
##João Lopes
set.seed(1)
setwd("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Ano\\Proteomics")
##Load a R package
library(limma)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(seqinr)
library(data.table)
library(tidyverse)
library(tximport)
library(readxl)
library(tidyr)
library(multiGSEA)
windowsFonts(Calibri = windowsFont("Arial"))

########################################Folders########################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create folders
dir.create("Files")
dir.create("Files/Tables")
dir.create("Plots")

dir.create("Plots/Volcano plots")
dir.create("Plots/Heatmap")
dir.create("Plots/Dotplot")


##Directories
pathData <- "Data/"
pathFiles <- "Files/"
pathTables <- "Files/Tables/"

pathVolcano <- "Plots/Volcano plots/"
pathHeatmap <- "Plots/Heatmap/"
pathDotplot <- "Plots/Dotplot/"

pathBothOmicsPathways <- "C:/Users/jdpl2/OneDrive/Ambiente de Trabalho/Mestrado/2º Ano/Both Omics/"
#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
########################################Folders########################################



#######################################Functions#######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

relativeExpression<-function(treatment1, treatment2, Design, fit){
  contrast<-paste(treatment1, "-", treatment2, sep = "")
  treatment1vtreatment2<-makeContrasts(contrasts = contrast, levels = Design)
  print(treatment1vtreatment2)
  lrt<-glmLRT(fit, contrast = treatment1vtreatment2)
  print(topTags(lrt))
  return(lrt)
}

volcanoPlot<-function(pathFiles, colour, legend, title, contrast, i){
  load(paste(pathFiles, "lrt.RData", sep = ""))
  print(colnames(DEP))
  signif<--log10(DEP[,paste("Pvalue-", contrast, sep = "")])
  plot(DEP[,paste("logFC-", contrast, sep = "")], signif, pch = ".", xlab = expression(log[2] * FC), ylab = expression(log[10] * FDR), main = NULL, 
       cex.lab = 1.5,
       panel.first={
         points(0, 0, pch = 16, cex = 1e6, col = "grey95")
         grid(col = "white", lty = 1)
       }
  )
  points(DEP[which(DEP[contrast] == 1), paste("logFC-", contrast, sep = "")], -log10(DEP[which(DEP[contrast] == 1), paste("Pvalue-", contrast, sep = "")]), pch = 21, cex = 1.5, col = "black",bg = colour[1,])
  points(DEP[which(DEP[contrast] == -1), paste("logFC-", contrast, sep = "")], -log10(DEP[which(DEP[contrast] == -1), paste("Pvalue-", contrast, sep = "")]), pch = 21, cex = 1.5, col = "black",bg = colour[2,])
  #legend("bottomleft", inset = c(-0.06,-0.5), xpd = TRUE, pch = 20, bty = "n",
         #col = c(colour[1,], colour[2,]), 
         #legend = c(legend[1,], legend[2,]))
  mtext(side = 3,
        LETTERS[i],
        at =-19.6+(-i+1)*3, 
        line = 2,
        cex = 1.5
  )
}

heatmapGraph<-function(heatData,colCutoff, rowCutoff){
  heatData<-as.matrix(depTable[,c(3,7,11,15,19)])
  rownames(heatData)<-depTable$Description
  head(heatData)
  clustFunction <- function(x) hclust(x, method="complete")
  distFunction <- function(x) dist(x,method="euclidean")
  #clusterint (sidebars)
  cCol<-colorRampPalette(brewer.pal(9, "Set1"))
  colFit<-clustFunction(distFunction(t(heatData)))
  cClusters<-cutree(colFit,h=colCutoff)
  cHeight<-length(unique(as.vector(cClusters)))
  colHeight = cCol(cHeight)
  cRow<-colorRampPalette(brewer.pal(9,"Set1"))
  rowFit<-clustFunction(distFunction(heatData))
  rClusters<-cutree(rowFit,h=rowCutoff)
  rHeight<-length(unique(as.vector(rClusters)));
  rowHeight = cRow(rHeight)
  xDimension<-25
  yDimension<-25
  windows(xDimension, yDimension)
  par(family="Arial")
  #Heatmap
  heatColour<-colorRampPalette(c("#fff04a", "black", "#c23e40"))(n = 100)
  p <- heatmap.2(heatData,
            hclust=clustFunction, 
            distfun=distFunction, 
            ColSideColors=colHeight[cClusters],
            RowSideColors=rowHeight[rClusters],
            density.info="none",
            col=heatColour, 
            trace="none", 
            scale="row", 
            cexCol = 0.8,
            celRow = 1.5,
            lhei = c(1,5),
            lwid = c(1.5,5), 
            offsetRow = 0, 
            offsetCol = 0,
            srtCol = 45, 
            key.title = NA, 
            margins = c(8,21.75),
            labRow = rownames(heatData))
  recordPlot()
  dev.print(tiff, filename = paste(pathHeatmap, "Heatmap", ".tif", sep = ""), height = 15, width = 15, units = "cm", res = 600)
}

dotplot <- function (data,path) {
  o <- ggplot(data, aes(x = contrast, y = pathway)) +
    geom_point(aes(size = 1-padj, color = NES, fill = NES), shape = 21, stroke = 0.5) +  
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  
    scale_size_continuous(range = c(3, 15)) +
    labs(x = "Contrast", y = "Pathway", 
         size = "1-p.adjust", 
         color = "NES",
         fill = "NES") +  
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),  
          legend.position = "bottom",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          panel.grid.major = element_line(linewidth = 1.2, color = "grey90"),  
          panel.grid.minor = element_line(linewidth = 0.8, color = "grey90"),  
          panel.background = element_rect(fill = "grey", color = NA))  
  ggsave(filename = paste(path, "Prot_Dotplot.svg", sep = ""), plot = o, device = "svg", width = 15, height = 15)
  ggsave(filename = paste(path, "Prot_Dotplot.tiff", sep = ""), plot = o, device = "tiff", width = 15, height = 15)
}

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Functions#######################################



##################################Objects/Descriptions#################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

titleCTL_10vMHW2_10<-expression("CTL day 10 vs. MHW2 day 10")
titleCTL_25vMHW1_25<-expression("CTL day 25 vs. MHW1 day 25")
titleCTL_25vMHW2_25<-expression("CTL day 25 vs. MHW2 day 25")
titleMHW1_25vMHW2_25<-expression("MHW1 day 25 vs. MHW2 day 25")
titleMHW2_10vMHW2_25<-expression("MHW2 day 10 vs. MHW2 day 25")

legend<-data.frame(CTL_10vMHW2_10 = c("Overexpressed in the control group in day 10", "Overexpressed in the MHW2 group in day 10"),
                   CTL_25vMHW1_25 = c("Overexpressed in the control group in day 25", "Overexpressed in the MHW1 group in day 25"),
                   CTL_25vMHW2_25 = c("Overexpressed in the control group in day 25","Overexpressed in MHW2 group in day 25"),
                   MHW1_25vMHW2_25 = c("Overexpressed in the MHW1 group in day 25","Overexpressed in MHW2 group in day 25"),
                   MHW2_10vMHW2_25 = c("Overexpressed in the MHW2 group in day 10","Overexpressed in MHW2 group in day 25"))

contrasts <- c("CTL_10vMHW2_10","CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25","MHW2_10vMHW2_25")

colour<-data.frame(CTL_10vMHW2_10 = c("#0000FF", "#FFA500"), 
                   CTL_25vMHW1_25 = c("#800080", "#FFFF00"),
                   CTL_25vMHW2_25 = c("#FF0000", "#00FFFF"),
                   MHW1_25vMHW2_25 = c("#8B4513", "#ADFF2F"),
                   MHW2_10vMHW2_25 = c("#009900", "#800000"))

Levels <- c(
  "CTL_10",
  "MHW2_10",
  "CTL_10",
  "MHW2_10",
  "CTL_10",
  "MHW2_10",
  "CTL_25",
  "MHW2_10",
  "CTL_25",
  "MHW1_25",
  "CTL_25",
  "MHW1_25",
  "MHW2_25",
  "MHW1_25",
  "MHW2_25",
  "MHW1_25",
  "MHW2_25",
  "MHW2_25"
)

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
##################################Objects/Descriptions#################################



#######################################Files/Data######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##File with all hits
file <- read_xlsx(paste(pathData,"CPMSF_GPPF-CM-47_EACMPC1-18_pomatochistus_microps_TMT18_HiRIEF_results.xlsx",sep = ""))
####################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Files/Data######################################



######Reading the file and creating a matrix with only the normalized abundances#######
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

#Create dataframe with UniProt Accession code, Protein name and normalized abundance of each protein
abundance_matrix <- file[,c(6,7,26:43)]
####################################################################################################

#Remove NAs
abundance_matrix <- na.omit(abundance_matrix)
###########

#Extract the protein name from the description
for (i in 1:nrow(abundance_matrix)) {
  result <- unlist(strsplit(abundance_matrix$Description[i], "OS="))
  previous_part <- result[1]
  abundance_matrix$Description[i] <- previous_part
}
##############################################

#Column names
colnames(abundance_matrix) <- c("Accession", 
                                "Protein",
                                "Abundances Normalized T1CTLPm1 MD10(M) - 126",
                                "Abundances Normalized T3MHW2Pm1 MD10(M) - 127_N",
                                "Abundances Normalized T10CTLPm1 MD10(M) - 127_C",
                                "Abundances Normalized T5MHW2Pm1 MD10(M) - 128_N",
                                "Abundances Normalized T15CTLPm1 MD10(M) - 128_C",
                                "Abundances Normalized T7MHW2Pm1 MD10(M) - 129_N",
                                "Abundances Normalized T1CTLPm1 MD25(M) - 129_C",
                                "Abundances Normalized T13MHW2Pm1 MD10(M) - 130_N",
                                "Abundances Normalized T4CTLPm1 MD25(M) - 130_C",
                                "Abundances Normalized T2MHW1Pm1 MD25(M) - 131_N",
                                "Abundances Normalized T10CTLPm1 MD25(M) - 131_C",
                                "Abundances Normalized T6MHW1Pm1 MD25(M) - 132_N",
                                "Abundances Normalized T3MHW2Pm1 MD25(M) - 132_C",
                                "Abundances Normalized T9MHW1Pm1 MD25(M) - 133_N",
                                "Abundances Normalized T5MHW2Pm1 MD25(M) - 133_C",
                                "Abundances Normalized T12MHW1Pm1 MD25(M) - 134_N",
                                "Abundances Normalized T7MHW2Pm1 MD25(M) - 134_C",
                                "Abundances Normalized T13MHW2Pm1 MD25(M) - 135")
#############

#Nº of duplicates
sum(duplicated(abundance_matrix$Accession)) 
#################

#Duplicate with higher mean
w <- c()
for (i in c(1:length(unique(abundance_matrix$Accession)))){
  dups <- which(abundance_matrix$Accession %in% unique(abundance_matrix$Accession)[i])
  v <- c()
  for (j in dups){
    v <- c(v,mean(as.numeric(abundance_matrix[-c(1,2)][j,])))
  }
  w = c(w,dups[which(v == max(v))[1]])
}
abundance_matrix <- abundance_matrix[w,]
####################

#Save table with UniProt Accession code, Protein Name and Normalized Abundance for each protein
write.table(abundance_matrix, paste(pathFiles, "Full.csv", sep = ""), sep = ";", col.names = NA)
###############################################################################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
######Reading the file and creating a matrix with only the normalized abundances#######



#######################################Statistics######################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Create a DGEList object from a table of counts
Data <-
  DGEList(counts = abundance_matrix[3:20],
          group = Levels,
          ids = rownames(abundance_matrix))
################################################

##Create the generalized linear models for comparing the expression levels from the different treatments
Design <- model.matrix( ~ 0 + Levels, data = Data$samples)
colnames(Design) <- levels(Data$samples$group)
Data <- estimateDisp(Data, Design)
fit <- glmFit(Data, Design)
########################################################################################################

##Expression levels between the treatments
CTL_10vMHW2_10<-relativeExpression(colnames(Design)[1],colnames(Design)[4],Design,fit)
MHW2_10vMHW2_25<-relativeExpression(colnames(Design)[4],colnames(Design)[5],Design,fit)
CTL_25vMHW1_25<-relativeExpression(colnames(Design)[2],colnames(Design)[3],Design,fit)
CTL_25vMHW2_25<-relativeExpression(colnames(Design)[2],colnames(Design)[5],Design,fit)
MHW1_25vMHW2_25<-relativeExpression(colnames(Design)[3],colnames(Design)[5],Design,fit)
##########################################

##Merge the logFC,logCPM,LR and Pvalue with the proteins unique ids
Results<-cbind(CTL_10vMHW2_10$table, MHW2_10vMHW2_25$table, CTL_25vMHW1_25$table, CTL_25vMHW2_25$table, MHW1_25vMHW2_25$table)
Results<-cbind(abundance_matrix$Accession,abundance_matrix$Protein,Results)
colnames(Results)<-c(
  "Accession",
  "Description",
  "logFC-CTL_10vMHW2_10",
  "logCPM-CTL_10vMHW2_10",
  "LR-CTL_10vMHW2_10",
  "Pvalue-CTL_10vMHW2_10",
  "logFC-MHW2_10vMHW2_25",
  "logCPM-MHW2_10vMHW2_25",
  "LR-MHW2_10vMHW2_25",
  "Pvalue-MHW2_10vMHW2_25",
  "logFC-CTL_25vMHW1_25",
  "logCPM-CTL_25vMHW1_25",
  "LR-CTL_25vMHW1_25",
  "Pvalue-CTL_25vMHW1_25",
  "logFC-CTL_25vMHW2_25",
  "logCPM-CTL_25vMHW2_25",
  "LR-CTL_25vMHW2_25",
  "Pvalue-CTL_25vMHW2_25",
  "logFC-MHW1_25vMHW2_25",
  "logCPM-MHW1_25vMHW2_25",
  "LR-MHW1_25vMHW2_25",
  "Pvalue-MHW1_25vMHW2_25"
)
head(Results)
write.table(Results, paste(pathFiles, "Results.csv", sep = ""), sep = ";", col.names = NA)
###############################################

##Calculating the nº of DEPs in general 
expressionTable<-decideTests(Results[,grepl("Pvalue-", colnames(Results))],coefficients = Results[,grepl("logFC", colnames(Results))], adjust.method = "fdr", lfc = 1.5)
expressionTable<-as.data.frame(expressionTable)
colnames(expressionTable)<-c("CTL_10vMHW2_10", "MHW2_10vMHW2_25", "CTL_25vMHW1_25","CTL_25vMHW2_25","MHW1_25vMHW2_25")
DEP<-cbind(Results,expressionTable)
depTable<-DEP[which(abs(DEP$CTL_10vMHW2_10) == 1 | abs(DEP$MHW2_10vMHW2_25) == 1 | abs(DEP$CTL_25vMHW1_25) == 1 | abs(DEP$CTL_25vMHW2_25) == 1 | abs(DEP$MHW1_25vMHW2_25) == 1),]
nrow(depTable) #Nº of DEPs
write.table(depTable, paste(pathBothOmicsPathways, "DEP.csv", sep = ""), sep = ";", col.names = NA)
save(CTL_10vMHW2_10, MHW2_10vMHW2_25, CTL_25vMHW1_25, CTL_25vMHW2_25, MHW1_25vMHW2_25, file = "Files/lrt.RData")
#######################################

##Nº of DEPs per treatment
DEP_CTL_10vMHW2_10<-sum(depTable$CTL_10vMHW2_10 != 0);DEP_CTL_10vMHW2_10
DEP_MHW2_10vMHW2_25<-sum(depTable$MHW2_10vMHW2_25 != 0);DEP_MHW2_10vMHW2_25
DEP_CTL_25vMHW2_25<-sum(depTable$CTL_25vMHW2_25 != 0);DEP_CTL_25vMHW2_25
DEP_CTL_25vMHW1_25<-sum(depTable$CTL_25vMHW1_25 != 0);DEP_CTL_25vMHW1_25
DEP_MHW1_25vMHW2_25<-sum(depTable$MHW1_25vMHW2_25 != 0);DEP_MHW1_25vMHW2_25
##########################

##Nº of underexpressed proteins per treatment 
uedep_CTL_10vMHW2_10<-sum(depTable$CTL_10vMHW2_10 == -1);uedep_CTL_10vMHW2_10
uedep_MHW2_10vMHW2_25<-sum(depTable$MHW2_10vMHW2_25 == -1);uedep_MHW2_10vMHW2_25
uedep_CTL_25vMHW2_25<-sum(depTable$CTL_25vMHW2_25 == -1);uedep_CTL_25vMHW2_25
uedep_CTL_25vMHW1_25<-sum(depTable$CTL_25vMHW1_25 == -1);uedep_CTL_25vMHW1_25
uedep_MHW1_25vMHW2_25<-sum(depTable$MHW1_25vMHW2_25 == -1);uedep_MHW1_25vMHW2_25
#############################################

##Nº of overexpressed proteins per treatment 
oedep_CTL_10vMHW2_10<-sum(depTable$CTL_10vMHW2_10 == 1);oedep_CTL_10vMHW2_10
oedep_MHW2_10vMHW2_25<-sum(depTable$MHW2_10vMHW2_25 == 1);oedep_MHW2_10vMHW2_25
oedep_CTL_25vMHW2_25<-sum(depTable$CTL_25vMHW2_25 == 1);oedep_CTL_25vMHW2_25
oedep_CTL_25vMHW1_25<-sum(depTable$CTL_25vMHW1_25 == 1);oedep_CTL_25vMHW1_25
oedep_MHW1_25vMHW2_25<-sum(depTable$MHW1_25vMHW2_25 == 1);oedep_MHW1_25vMHW2_25
############################################

##Table with total DEP, underexpressed proteins and overexpressed proteins per treatment
dep_per_treatment <- data.frame(Treatments = contrasts, 
                                #No_of_Proteins = c(DEP_CTL_10vMHW2_10,DEP_CTL_25vMHW1_25,DEP_CTL_25vMHW2_25,DEP_MHW1_25vMHW2_25,DEP_MHW2_10vMHW2_25),
                                No_of_underexpressed_Proteins = c(uedep_CTL_10vMHW2_10,uedep_CTL_25vMHW1_25,uedep_CTL_25vMHW2_25,uedep_MHW1_25vMHW2_25,uedep_MHW2_10vMHW2_25),
                                No_of_overexpressed_Proteins = c(oedep_CTL_10vMHW2_10,oedep_CTL_25vMHW1_25,oedep_CTL_25vMHW2_25,oedep_MHW1_25vMHW2_25,oedep_MHW2_10vMHW2_25))
colnames(dep_per_treatment) <- c("Treatments","Nº of underexpressed proteins","Nº of overexpressed proteins")
write.table(dep_per_treatment, paste(pathTables, "dep_per_treatment.csv",sep=""), sep=";", col.names=TRUE, row.names = FALSE)
########################################################################################

##Create tables with accession, logFC, Pvalue and Description for all proteins
depTable_up <- DEP[which(DEP$CTL_10vMHW2_10 == 1 | DEP$MHW2_10vMHW2_25 == 1 | DEP$CTL_25vMHW1_25 == 1 | DEP$CTL_25vMHW2_25 == 1 | DEP$MHW1_25vMHW2_25 == 1),]
##############################################################################

##Create tables with accession, logFC, Pvalue and Description for each treatment
depTable_CTL_10vMHW2_10 <- DEP[which(abs(DEP$CTL_10vMHW2_10) == 1),]
depTable_CTL_10vMHW2_10 <- depTable_CTL_10vMHW2_10[, c(1,2,3,6)]
depTable_CTL_10vMHW2_10$padj <- p.adjust(depTable_CTL_10vMHW2_10$`Pvalue-CTL_10vMHW2_10`, method = "fdr")
write.table(depTable_CTL_10vMHW2_10, paste(pathTables, "dep_CTL_10VMHW2_10.csv", sep = ""), sep = ";", col.names = TRUE, row.names = FALSE)

depTable_MHW2_10vMHW2_25 <- DEP[which(abs(DEP$MHW2_10vMHW2_25) == 1),]
depTable_MHW2_10vMHW2_25 <- depTable_MHW2_10vMHW2_25[, c(1,2,7,10)]
depTable_MHW2_10vMHW2_25$padj <- p.adjust(depTable_MHW2_10vMHW2_25$`Pvalue-MHW2_10vMHW2_25`, method = "fdr")
write.table(depTable_MHW2_10vMHW2_25, paste(pathTables, "dep_MHW2_10vMHW2_25.csv", sep = ""), sep = ";", col.names = TRUE, row.names = FALSE)

depTable_CTL_25vMHW2_25 <- DEP[which(abs(DEP$CTL_25vMHW2_25) == 1),]
depTable_CTL_25vMHW2_25 <- depTable_CTL_25vMHW2_25[, c(1,2,15,18)]
depTable_CTL_25vMHW2_25$padj <- p.adjust(depTable_CTL_25vMHW2_25$`Pvalue-CTL_25vMHW2_25`, method = "fdr")
write.table(depTable_CTL_25vMHW2_25, paste(pathTables, "dep_CTL_25vMHW2_25.csv", sep = ""), sep = ";", col.names = TRUE, row.names = FALSE)

depTable_CTL_25vMHW1_25 <- DEP[which(abs(DEP$CTL_25vMHW1_25) == 1),]
depTable_CTL_25vMHW1_25 <- depTable_CTL_25vMHW1_25[, c(1,2,11,14)]
depTable_CTL_25vMHW1_25$padj <- p.adjust(depTable_CTL_25vMHW1_25$`Pvalue-CTL_25vMHW1_25`, method = "fdr")
write.table(depTable_CTL_25vMHW1_25, paste(pathTables, "dep_CTL_25vMHW1_25.csv", sep = ""), sep = ";", col.names = TRUE, row.names = FALSE)

depTable_MHW1_25vMHW2_25 <- DEP[which(abs(DEP$MHW1_25vMHW2_25) == 1),]
depTable_MHW1_25vMHW2_25 <- depTable_MHW1_25vMHW2_25[, c(1,2,19,22)]
depTable_MHW1_25vMHW2_25$padj <- p.adjust(depTable_MHW1_25vMHW2_25$`Pvalue-MHW1_25vMHW2_25`, method = "fdr")
write.table(depTable_MHW1_25vMHW2_25, paste(pathTables, "dep_MHW1_25vMHW2_25.csv", sep = ""), sep = ";", col.names = TRUE, row.names = FALSE)
################################################################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#######################################Statistics######################################



####################################Pathway analysis###################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

DEP_pathway <- DEP[c(1,3:22)]
write.table(DEP_pathway, paste(pathBothOmicsPathways, "DEP_pathway.csv", sep = ""),sep = ";",col.names = TRUE, row.names = FALSE)

##Create data frames with Accession, logFC and Pvalue for the pathway analysis (all proteins)
allp_CTL_10vMHW2_10<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC-CTL_10vMHW2_10`,DEP_pathway$`Pvalue-CTL_10vMHW2_10`)
names(allp_CTL_10vMHW2_10) <- c("Accession","logFC","Pvalue")

allp_CTL_25vMHW1_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC-CTL_25vMHW1_25`,DEP_pathway$`Pvalue-CTL_25vMHW1_25`)
names(allp_CTL_25vMHW1_25) <- c("Accession","logFC","Pvalue")

allp_CTL_25vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC-CTL_25vMHW2_25`,DEP_pathway$`Pvalue-CTL_25vMHW2_25`)
names(allp_CTL_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allp_MHW1_25vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC-MHW1_25vMHW2_25`,DEP_pathway$`Pvalue-MHW1_25vMHW2_25`)
names(allp_MHW1_25vMHW2_25) <- c("Accession","logFC","Pvalue")

allp_MHW2_10vMHW2_25<- data.frame(DEP_pathway$Accession,DEP_pathway$`logFC-MHW2_10vMHW2_25`,DEP_pathway$`Pvalue-MHW2_10vMHW2_25`)
names(allp_MHW2_10vMHW2_25) <- c("Accession","logFC","Pvalue")
############################################################################

##Create a data structure
layers <- c("proteome")
odataCTL_10vMHW2_10 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW1_25 <- initOmicsDataStructure(layer=layers)
odataCTL_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW1_25vMHW2_25 <- initOmicsDataStructure(layer=layers)
odataMHW2_10vMHW2_25 <- initOmicsDataStructure(layer=layers)
#########################

##Add proteome layer
odataCTL_10vMHW2_10$proteome <- rankFeatures(allp_CTL_10vMHW2_10$logFC,allp_CTL_10vMHW2_10$Pvalue)
names(odataCTL_10vMHW2_10$proteome) <- allp_CTL_10vMHW2_10$Accession
odataCTL_10vMHW2_10$proteome <- sort(odataCTL_10vMHW2_10$proteome)
head(odataCTL_10vMHW2_10$proteome)

odataCTL_25vMHW1_25$proteome <- rankFeatures(allp_CTL_25vMHW1_25$logFC,allp_CTL_25vMHW1_25$Pvalue)
names(odataCTL_25vMHW1_25$proteome) <- allp_CTL_25vMHW1_25$Accession
odataCTL_25vMHW1_25$proteome <- sort(odataCTL_25vMHW1_25$proteome)
head(odataCTL_25vMHW1_25$proteome)

odataCTL_25vMHW2_25$proteome <- rankFeatures(allp_CTL_25vMHW2_25$logFC,allp_CTL_25vMHW2_25$Pvalue)
names(odataCTL_25vMHW2_25$proteome) <- allp_CTL_25vMHW2_25$Accession
odataCTL_25vMHW2_25$proteome <- sort(odataCTL_25vMHW2_25$proteome)
head(odataCTL_25vMHW2_25$proteome)

odataMHW1_25vMHW2_25$proteome <- rankFeatures(allp_MHW1_25vMHW2_25$logFC,allp_MHW1_25vMHW2_25$Pvalue)
names(odataMHW1_25vMHW2_25$proteome) <- allp_MHW1_25vMHW2_25$Accession
odataMHW1_25vMHW2_25$proteome <- sort(odataMHW1_25vMHW2_25$proteome)
head(odataMHW1_25vMHW2_25$proteome)

odataMHW2_10vMHW2_25$proteome <- rankFeatures(allp_MHW2_10vMHW2_25$logFC,allp_MHW2_10vMHW2_25$Pvalue)
names(odataMHW2_10vMHW2_25$proteome) <- allp_MHW2_10vMHW2_25$Accession
odataMHW2_10vMHW2_25$proteome <- sort(odataMHW2_10vMHW2_25$proteome)
head(odataMHW2_10vMHW2_25$proteome)
####################

##Select the databases we want to query and download pathway definitions
databases <- c("kegg")
pathways <- getMultiOmicsFeatures(dbs = databases, layer = layers,
                                  returnProteome = "UNIPROT",
                                  organism = "drerio",
                                  useLocal =  FALSE)
pathways_short <- lapply(names(pathways), function(name) {
  head(pathways[[name]], 2)
})
names(pathways_short) <- names(pathways)
pathways_short
pathways$proteome[8]
########################################################################

##Run the pathway enrichment
enrichment_scoresCTL_10vMHW2_10 <- multiGSEA(pathways,odataCTL_10vMHW2_10)
Tenrichment_scoresCTL_10vMHW2_10 <- as.data.frame(enrichment_scoresCTL_10vMHW2_10$proteome)
Tenrichment_scoresCTL_10vMHW2_10$leadingEdge <- sapply(Tenrichment_scoresCTL_10vMHW2_10$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_10vMHW2_10 <- Tenrichment_scoresCTL_10vMHW2_10[order(Tenrichment_scoresCTL_10vMHW2_10$padj),]
Tenrichment_scoresCTL_10vMHW2_10$pathway <- sub("^\\(KEGG\\) ", "", Tenrichment_scoresCTL_10vMHW2_10$pathway)
top_10_Tenrichment_scoresCTL_10vMHW2_10 <- head(Tenrichment_scoresCTL_10vMHW2_10,10)
top_10_Tenrichment_scoresCTL_10vMHW2_10$contrast <- rep("CTL_10vMHW2_10")
write.table(top_10_Tenrichment_scoresCTL_10vMHW2_10,paste(pathTables,"top10_enrichment_scoresCTL_10vMHW2_10.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW1_25 <- multiGSEA(pathways,odataCTL_25vMHW1_25)
Tenrichment_scoresCTL_25vMHW1_25 <- as.data.frame(enrichment_scoresCTL_25vMHW1_25$proteome)
Tenrichment_scoresCTL_25vMHW1_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW1_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW1_25 <- Tenrichment_scoresCTL_25vMHW1_25[order(Tenrichment_scoresCTL_25vMHW1_25$padj),]
Tenrichment_scoresCTL_25vMHW1_25$pathway <- sub("^\\(KEGG\\) ", "", Tenrichment_scoresCTL_25vMHW1_25$pathway)
top_10_Tenrichment_scoresCTL_25vMHW1_25 <- head(Tenrichment_scoresCTL_25vMHW1_25,10)
top_10_Tenrichment_scoresCTL_25vMHW1_25$contrast <- rep("CTL_25vMHW1_25")
write.table(top_10_Tenrichment_scoresCTL_25vMHW1_25,paste(pathTables,"top10_enrichment_scoresCTL_25vMHW1_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresCTL_25vMHW2_25 <- multiGSEA(pathways,odataCTL_25vMHW2_25)
Tenrichment_scoresCTL_25vMHW2_25 <- as.data.frame(enrichment_scoresCTL_25vMHW2_25$proteome)
Tenrichment_scoresCTL_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresCTL_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresCTL_25vMHW2_25 <- Tenrichment_scoresCTL_25vMHW2_25[order(Tenrichment_scoresCTL_25vMHW2_25$padj),]
Tenrichment_scoresCTL_25vMHW2_25$pathway <- sub("^\\(KEGG\\) ", "", Tenrichment_scoresCTL_25vMHW2_25$pathway)
top_10_Tenrichment_scoresCTL_25vMHW2_25 <- head(Tenrichment_scoresCTL_25vMHW2_25,10)
top_10_Tenrichment_scoresCTL_25vMHW2_25$contrast <- rep("CTL_25vMHW2_25")
write.table(top_10_Tenrichment_scoresCTL_25vMHW2_25,paste(pathTables,"top10_enrichment_scoresCTL_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW1_25vMHW2_25 <- multiGSEA(pathways,odataMHW1_25vMHW2_25)
Tenrichment_scoresMHW1_25vMHW2_25 <- as.data.frame(enrichment_scoresMHW1_25vMHW2_25$proteome)
Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW1_25vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW1_25vMHW2_25 <- Tenrichment_scoresMHW1_25vMHW2_25[order(Tenrichment_scoresMHW1_25vMHW2_25$padj),]
Tenrichment_scoresMHW1_25vMHW2_25$pathway <- sub("^\\(KEGG\\) ", "", Tenrichment_scoresMHW1_25vMHW2_25$pathway)
top_10_Tenrichment_scoresMHW1_25vMHW2_25 <- head(Tenrichment_scoresMHW1_25vMHW2_25,10)
top_10_Tenrichment_scoresMHW1_25vMHW2_25$contrast <- rep("MHW1_25vMHW2_25")
write.table(top_10_Tenrichment_scoresMHW1_25vMHW2_25,paste(pathTables,"top10_enrichment_scoresMHW1_25vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)

enrichment_scoresMHW2_10vMHW2_25 <- multiGSEA(pathways,odataMHW2_10vMHW2_25)
Tenrichment_scoresMHW2_10vMHW2_25 <- as.data.frame(enrichment_scoresMHW2_10vMHW2_25$proteome)
Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge <- sapply(Tenrichment_scoresMHW2_10vMHW2_25$leadingEdge, function(x) paste(x, collapse = ";"))
Tenrichment_scoresMHW2_10vMHW2_25 <- Tenrichment_scoresMHW2_10vMHW2_25[order(Tenrichment_scoresMHW2_10vMHW2_25$padj),]
Tenrichment_scoresMHW2_10vMHW2_25$pathway <- sub("^\\(KEGG\\) ", "", Tenrichment_scoresMHW2_10vMHW2_25$pathway)
top_10_Tenrichment_scoresMHW2_10vMHW2_25 <- head(Tenrichment_scoresMHW2_10vMHW2_25,10)
top_10_Tenrichment_scoresMHW2_10vMHW2_25$contrast <- rep("MHW2_10vMHW2_25")
write.table(top_10_Tenrichment_scoresMHW2_10vMHW2_25,paste(pathTables,"top10_enrichment_scoresMHW2_10vMHW2_25.csv",sep=""),sep=";",row.names = FALSE)
############################

##Create dataset with all enrichment scores
combined_df <- rbind(top_10_Tenrichment_scoresCTL_10vMHW2_10, top_10_Tenrichment_scoresCTL_25vMHW1_25, top_10_Tenrichment_scoresCTL_25vMHW2_25, top_10_Tenrichment_scoresMHW1_25vMHW2_25, top_10_Tenrichment_scoresMHW2_10vMHW2_25)
###########################################

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
####################################Pathway analysis###################################



##########################################Plots########################################
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV#

##Volcano Plot
par(mfrow = c(1,1), mar = c(7,7,4,4), family = "Arial")
for(i in 1:length(contrasts)){
  volcanoPlot(pathFiles, colour[contrasts[i]], legend[contrasts[i]], get(paste("title", contrasts[i], sep = "")), contrasts[i], i)
  dev.print(tiff, filename = paste(pathVolcano, "VolcanoPlot-", contrasts[i], ".tif", sep = ""), height = 15, width = 15, units = "cm", res = 600)
  dev.copy(svg, filename = paste(pathVolcano, "VolcanoPlot-",contrasts[i], ".svg", sep = ""), height = 15, width = 15)
  dev.off()
}
##############

##Heatmap
colCutoff = 10
rowCutoff = 4
heatmapGraph(heatData,colCutoff, rowCutoff)
#########

#Dotplot
dotplot(combined_df,pathDotplot)
########

#ΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛΛ#
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#
##########################################Plots########################################