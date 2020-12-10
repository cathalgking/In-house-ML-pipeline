###############################################################################################################
# This script will:
# 1. Take as input a master matrix and a design marix and Normalize the master matrix by analysing a PDF report.
# 2. Take the normalized matrix of choice and calculate statistics (log2FC, p-value and PDF report).
# 3. Make volcano plots from any 2 chosen groups in the normalised dataset.
# 4. Collect significant genes from 'top_left' and 'top_right' areas of the plot. 
###############################################################################################################

### NormalyzerDE ###

#load relevant libraries
library(NormalyzerDE)
library(dplyr)
library(matrixStats)
library(EnhancedVolcano)

# outDir should contain the PATH to where your output report and normalised matrices will be 
outDir <- ('/home/cathal_king/Desktop')

# read in, write out
#df2 <- read.table(file = "/home/cathal_king/Desktop/CycLoess-normalized_round2_zeros.txt", header = TRUE)
#write.table(df2, file = "/home/cathal_king/Desktop/CycLoess-normalized_round2_zeros2.txt", row.names = FALSE, sep = '\t', quote = FALSE)

# The following line reads in the design matrix. In this example, my design matrix file (design_matrix1.csv) is in the folder "FYP_datasets" which is in my current working directory.
# Set your working directory with setwd() and check what it is with getwd()
designFp <- "/home/cathal_king/Downloads/design_matrix_round2.txt"

# Then read in your master dataset
dataFp <- "/home/cathal_king/Desktop/TBI_2/TBI2_out_dedup4.txt"

# Finally, call the NormalyzerDE command with the following line. You should then see it processing in the Console. 
# Add any other relevant Arguments to the commands by checking the package notes.
# jobName refers the folder name that will be created in the outDir.
NormalyzerDE::normalyzer(jobName="Cathals_TBI2run", designPath = designFp, dataPath = dataFp, outputDir=outDir, tinyRunThres = 10000000, noLogTransform = T, omitLowAbundSamples = T, sampleAbundThres = 1110)
# noLogTransform = T, omitLowAbundSamples = T, sampleAbundThres = 110, , skipAnalysis = T
## Notes ###
# Both your design and data matrices ideally should be tab seperated text files or Text CSV files. Comma seperated will not work.
# Calling the arguments designFp or dataFp in the console will not show a matrix. These arguments are designed to take file paths directly, rather than load data frames.
# Run "??NormalyzerDE" in the Console to explore some other useful arguments that the package provides and can be added to this script. 

# statistical report
# import a design matrix with exact numberings for the statistical analyses below or use the one above as long as they have the correct groupings for the volcano plots
# designFp <- "/home/cathal_king/Desktop/test_Design_matrix.txt"
# run the stats report and choose the groupings in comparison=c() . Choose the normalised matrix and input next line
normMatrixPath <- paste(outDir, "CycLoess-normalized_round2_zeros3.txt", sep="/")

# this next command will carry out statistical analyses (pvalue & logFC) on the normalised dataset
normalyzerDE("round2_stats_run38-36",
             comparisons=c("38-36"),
             designPath=designFp,
             dataPath=normMatrixPath,
             outputDir=outDir,
             condCol="group")



### EnhancedVolcano ###

### Volcano Plot

# read in house keeping genes list
hkg <- read.csv('/home/cathal_king/Downloads/HK_genes_only.txt')


# read in the dataset containing the extra statistics columns created from NormalyzerDE above
res <- read.table("/home/cathal_king/Desktop/round2_stats_run38-36/round2_stats_run38-36_stats.tsv", sep = '\t', header = T)

# inter quartile range (IQR). Find the IQR for the HKG only. First, remove all other genes in the dataset.

## remove rows that dont have hk genes.
res_only_hkg <- hkg %>% inner_join(res, by = 'GENE')

# either export that matrix to Perseus to calculate the IQR or do it here in R
#write.table(res_only_hkg, file = "res_hkg2_5-6.txt", row.names = FALSE, sep = ",")

# get the IQR of the column
# remove the p-value and fold change columns (the 2nd to the 5th columns)
res_only_hkg_just_data <- select(res_only_hkg, -c(1:5))


## IMPORTANT ## - change names of columns to keep in the next line of code
# only keep the columns of interest to calculate the IQR's. enter the first and last column names to keep
res_only_hkg_just_data_cols <- select(res_only_hkg_just_data, D24.G0.KIDNEY.1:D24.G0.LIVER.5) 
# get the IQR's of all the rows
mtx <- data.matrix(res_only_hkg_just_data_cols)
IQRS <- rowIQRs(mtx)
# get median of that list, and make a record of it
average = mean(IQRS)


# The EnhancedVolcano package expects a dataset containing p-value(adjusted or unadjusted), log2FC and the genes/identifiers of interest. 
# If you are importing data that already has been converted to -log10, for example perseus dataset (from a saved perseus session), they can be converted back to p-value with the folowing line of code. 
# res$pvalue <- 1/10^(res$X.Log.P.value.)
# res$X4.8_AdjPVal
# Generate the plot. x and y arguments refer to the columns in the dataset that contain the log fold change and p-values.
# Customise as needed according to the vignette found at https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
# x is the column containing FC and y is the column containing pvalue
EnhancedVolcano(res,
                lab = NA,
                x = 'X38.36_log2FoldChange',
                y = 'X38.36_PValue',
                xlim = c(-4, 4),
                ylim = c(0,35),
                title = "G.P.Finkielstain et al 2009(GSE38754) - Liver",
                pCutoff = 0.05,
                FCcutoff = 0.55,
                pointSize = 4,
                legendPosition = 'bottom',
                hlineType = 'longdash')


# collect genes from significant areas of the plot
# Important - tweak following commands (numbers after pval's and lFC's) to suit thresholds and groupings calculated above
top_right <- subset(res, X38.36_PValue < 0.05 & X38.36_log2FoldChange>0.55)
top_left <- subset(res, X38.36_PValue < 0.05 & X38.36_log2FoldChange<(0.55*-1))

write.csv(top_right$GENE, file = "top right 38-36 ECM Round 2", row.names = FALSE, quote = FALSE)
write.csv(top_left$GENE, file = "top left 38-36 ECM Round 2", row.names = FALSE, quote = FALSE)

  
# if genes in both significant areas are required then uncomment, change the FC number & run the following line
# all_sig_genes <- subset(res, X14.13_PValue < 0.05 & abs(X14.13_AdjPVal)>0.03)

