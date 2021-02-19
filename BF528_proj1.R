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
rma <- rma(all_data)

#normalize for QC
data = fitPLM(all_data, normalize = TRUE, background = TRUE)

#QC using NUSE, anticipating median approx. 1
NUSE_stats <- NUSE(data,type="stats")
hist(NUSE_stats[1, ])

#QC using RLE, anticipating median approx. 0
RLE_stats <- RLE(data,type="stats")
hist(RLE_stats[1, ])

#import metadata for ComBat
file <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv', header = TRUE)
batch <- file$normalizationcombatbatch
mod <- model.matrix(~normalizationcombatmod, data=file)

#pickup here!!!
write.csv(exprs(rma), '/projectnb/bf528/users/test/test_expression.csv')
expression = (exprs(rma))
cb = ComBat(expression, batch, mod)
write.csv(cb, '/projectnb/bf528/users/test/test_expression_cb.csv')

scaled = scale(expression)
pca = prcomp(scaled, scale = FALSE, center = FALSE)
pca_results = pca$rotation

PC1_imp = toString(summary(pca)$importance[2, 1])
PC2_imp = summary(pca)$importance[2, 2]

PC1_label <- paste(c("PC1: Importance =", PC1_imp), collapse = " ")
PC2_label <- paste(c("PC2: Importance =", PC2_imp), collapse = " ")

pca$rotation[,1:2]

plot(x = pca$rotation[,1:1], y = pca$rotation[,2:2] ,xlab = PC1_label, ylab = PC2_label)

summary(pca)$importance[2:2, ]

#just curious now...
sample_data = read.csv('/project/bf528/project_1/data/example_intensity_data.csv', sep=" ")
total_reads = dim(sample_data)[1]

exp_results <- read.csv("/project/bf528/project_1/data/differential_expression_results.csv")
matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(row.names(exp_results)),columns = ("SYMBOL"))
#merge exp_results w/ matches and drop duplicates based on 
exp_results$PROBE_ID <- rownames(exp_results)
exp_data <- merge(exp_results, matches, by.x = 4, by.y = 1, all = TRUE)


