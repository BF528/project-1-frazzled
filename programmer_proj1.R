#BiocManager::install(c('affy', 'affyPLM', 'sva', 'AnnotationDbi',  'hgu133plus2.db'))

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

#reading in files and including extra sample
CEL_files = ReadAffy(celfile.path ='/projectnb/bf528/users/frazzled/project_1/samples/CEL_files/' )
extra_sample = ReadAffy(celfile.path ='/projectnb/bf528/users/frazzled/project_1/samples/' )
all_data = merge(CEL_files, extra_sample)

#normalize data using RMA
rma_data <- rma(all_data)

#normalize for QC
normed_data = fitPLM(all_data, normalize = TRUE, background = TRUE)

#QC using NUSE, anticipating median approx. 1
NUSE_stats <- NUSE(normed_data,type="stats")
hist(NUSE_stats[1, ])

#QC using RLE, anticipating median approx. 0
RLE_stats <- RLE(normed_data,type="stats")
hist(RLE_stats[1, ])

#import metadata for ComBat
file <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv', header = TRUE)
batch <- file$normalizationcombatbatch
mod <- model.matrix(~normalizationcombatmod, data=file)

#adjusting expression wiht ComBat and exporting to CSV
expression = (exprs(rma_data))
#write.csv(expression, '/projectnb/bf528/users/frazzled/project_1.csv')
cb_adjusted_data = ComBat(expression, batch, mod)
write.csv(cb_adjusted_data, '/projectnb/bf528/users/frazzled/project_1.csv')

#scale data and create PCA
scaled = scale(expression)
pca = prcomp(scaled, scale = FALSE, center = FALSE)
pca_results = pca$rotation

#extract importance of each PC
PC1_imp = summary(pca)$importance[2, 1]
PC2_imp = summary(pca)$importance[2, 2]

#creating axis labels for graph
PC1_label <- paste(c("PC1: Importance =", PC1_imp), collapse = " ")
PC2_label <- paste(c("PC2: Importance =", PC2_imp), collapse = " ")

#plotting PC1 and PC2
plot(pca$rotation[,1:2] ,xlab = PC1_label, ylab = PC2_label)

#comparing my results to sample input for analyst, out of curiosity
#sample_data = read.csv('/project/bf528/project_1/data/example_intensity_data.csv', sep=" ")



