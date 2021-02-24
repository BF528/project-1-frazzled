#install hgu133plus2.db package 
BiocManager::install(c("hgu133plus2.db"))
#install dplyr and tidyverse
install.packages('dplyr')
install.packages('tidyverse')
#install GSEABase
BiocManager::install(c("GSEABase"))

#run libraries
library("hgu133plus2.db")
library(dplyr)
library(tidyverse)
library("GSEABase")

#read csv 
noise <- read.csv("/projectnb/bf528/users/frazzled/project_1/final_Noise_filter.csv")
results<- read.csv("/projectnb/bf528/users/frazzled/project_1/final_section5_result.csv") #5.4
n_results <-read.csv("/projectnb/bf528/users/frazzled/project_1/final_section5.6_result.csv")
#match PROBEIDs with SYMBOLS
matches <- AnnotationDbi::select(hgu133plus2.db, keys = as.character(n_results$X),columns = ("SYMBOL"))

#change the name of X in results to PROBEID 
results <- n_results %>% 
  rename(PROBEID = X)
#merge results and matches based on PROBEID
mapped_results <- left_join(results,matches)
#drop NAs
mapped_r <- drop_na(mapped_results)
mapped <- mapped_r %>% distinct()

# added column abso_t_value and calculated abs value of t then arranged by most expressed 
r_by_abst <- mutate(mapped,"abso_t_value"= abs(t.stat)) %>% 
  arrange(desc(abso_t_value))
#take top 1000 and top10 differentially expressed genes
diff_exp_top1000 <- r_by_abst[1:1000,]
diff_exp_top10 <- r_by_abst[1:10,]
not_diff_exp <- r_by_abst[1001:dim(r_by_abst)[1],"SYMBOL"]
down_reg_top10 <- not_diff_exp[10:1,]

#read in gmt files 
gscKEGG<-getGmt("c2.cp.kegg.v7.2.symbols.gmt")
gscGO<-getGmt("c5.go.v7.2.symbols.gmt")
gscHALL<-getGmt("h.all.v7.2.symbols.gmt")

#listing gene sets 
KEGG_GID <- diff_exp_top1000$SYMBOL
#looping through gscKEGG objects to pull out each sets geneids and saving geneids to v 
geneset_KEGG <- for( i in names(gscKEGG)) {
  v <- geneIds(gscKEGG[[i]])
  x <- intersect(KEGG_GID,v)  #intersect took the common genes between v and KEGG_GID x = genes in set that are differentially expressed
  in_geneset <- length(x)
  not_in_geneset <- length(KEGG_GID) - in_geneset #is differentially exp not in gene set 
  not_diff_exp_igs <- intersect(not_diff_exp,v) #not differentially exp in gene set 
  l <- length(not_diff_exp_igs)
  not_in_gs_not_diff_exp <- length(not_diff_exp) - l #not in gene set and not diff exp genes
}

GO_GID<- diff_exp_top1000$SYMBOL
#looping through gscGO to pull out geneids and saving gids
geneset_GO<- for( n in names(gscGO)) {
  j <- geneIds(gscGO[[n]])
  k <- intersect(GO_GID,j) #intersect getting common genes between each objects geneids in j and the genes in gscGO that are diff exp
  go_in_geneset <- length(k) #in geneset 
  go_notin_geneset <- length(GO_GID) - go_in_geneset #is diff exp not in gene set
  go_not_diff_exp_igs <- intersect(not_diff_exp,j) #not diff exp in gene set
  m <- length(go_not_diff_exp_igs)
  go_notigs_notdiffexp <- length(not_diff_exp) - m #not in gene set and not diff exp genes
}

HALL_GID<- diff_exp_top1000$SYMBOL
#looping through gscHALL to pull out geneids and saving 
geneset_HALL <- for (y in names(gscHALL)) {
  gid <-geneIds(gscHALL[[y]])
  h <- intersect(HALL_GID,gid) #getting common genes that are diff expressed
  hall_ings <- length(h)
  hall_not_ings <- length(HALL_GID) - hall_ings #diff exp genes not in gene set
  hall_notdiffexp_igs <- intersect(not_diff_exp,gid) #not diff expressed genes in geneset
  h_hall <- length(hall_notdiffexp_igs)
  hall_notigs_notdiffexp <- length(not_diff_exp) - h_hall #not diff expressed not in geneset 
}
#create contingency tables 
  #contingency table for KEGG 
fish_test_kegg <- fisher.test(matrix(c(6,994,68,46595),nrow = 2))
  #contingency table for GO
fish_test_go <- fisher.test(matrix(c(1,999,7,46656),nrow = 2))
  #contingency table for HALL 
fish_test_hall <- fisher.test(matrix(c(3,997,39,46624),nrow = 2))

#adj p 
adjp_kegg_ft <- p.adjust(fish_test_kegg$p.value,method = 'BH')
adjp_hall_ft <- p.adjust(fish_test_hall$p.value,method = 'BH')
adjp_go_ft <- p.adjust(fish_test_go$p.value,method = 'BH')

#significance 
kegg_sig <-adjp_kegg_ft<0.05
hall_sig <-adjp_hall_ft<0.05
go_sig <- adjp_go_ft<0.05

#selecting top 3 de genesets
