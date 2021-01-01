##This is the main R parsing script for the AF100 project 
##It performs basic data cleaning, visualization, and statistics.
#The main input files are the output of snpEff 
#and manually input clade designations with a timeline of n-days since isolation 

#load modules
require(data.table)
require(splitstackshape)
library(plyr)
library(dplyr)
library(tidyverse)
library(rlist)
library(gdata)

#set dir
setwd(<>)
options(stringsAsFactors = FALSE)

#load main dat files from snpEff
snpEff_all<-read.delim("Big_tree_clean.snpEff.matrix.tab", header = TRUE, sep = "\t", fill = TRUE, strip.white = TRUE)

###clean the data for easy processing###
#remove the brackets form ALT
snpEff_all$ALT<-gsub("\\[|\\]", "", snpEff_all$ALT)

#remove the "/" from all the the SNP calls 
snpEff_all[] <- lapply(snpEff_all, gsub, pattern='/', replacement="")

#get set size
dim(snpEff_all[,11:(ncol(snpEff_all) -1)])
#121 strains
#223369 variants

#remove intronic and synonamous variants 
no_intergenic_a<- snpEff_all[!snpEff_all$TYPE == "intergenic",]
no_intergenic_b<- no_intergenic_a[!no_intergenic_a$TYPE == "intron_variant",]
snpEff_no_intergenic<- no_intergenic_b[!no_intergenic_b$TYPE == "synonymous_variant",]

#check
#dim
dim(snpEff_all)
dim(snpEff_no_intergenic)

###set groups for all downstream analysis###
#Outgoups
OG_names<- c("DMC2_AF100.7B", 
             "DMC2_AF100.2B",
             "AF100.11B",
             "X1F1SW_F4")
#Clade 1
Clade_1_names<- c("DMC2_AF100.3B", 
                  "DMC2_AF100.4B", 
                  "DMC2_AF100.5B", 
                  "DMC2_AF100.6B", 
                  "AF100.10B",
                  "AF100.12_7")

#Clade 1 with ancestor
Clade_1_w_ancestor_names<- c("DMC2_AF100.1_3",
                             "DMC2_AF100.3B", 
                             "DMC2_AF100.4B", 
                             "DMC2_AF100.5B", 
                             "DMC2_AF100.6B", 
                             "AF100.10B",
                             "AF100.12_7")
#Clade 2
#note - there is no clear Clade 2 ancestor 
Clade_2_names<- c("DMC2_AF100.8B", 
                  "DMC2_AF100.9B", 
                  "AF100.10_5", 
                  "AF100.11_3", 
                  "DMC2_AF100.12_9")


#subset the clades from the main df
OG<- snpEff_no_intergenic[, OG_names]
Clade_1<- snpEff_no_intergenic[, Clade_1_names]
Clade_1_w_ancestor<- snpEff_no_intergenic[, Clade_1_w_ancestor_names]
Clade_2<- snpEff_no_intergenic[, Clade_2_names]

#subset the main df from the clades
All_but_clades1<- snpEff_no_intergenic[,12:ncol(snpEff_no_intergenic) -1]
All_but_clades2<- All_but_clades1[!names(All_but_clades1) %in% Clade_1_names]
All_but_clades<- All_but_clades2[!names(All_but_clades2) %in% Clade_2_names]

#check - should be 123 - 11 = 112
dim(snpEff_no_intergenic[,11:(ncol(snpEff_no_intergenic) -1)])
dim(All_but_clades)

#calculate new ALT for each group (calculates the alternative allele for the clade: i.e. if the alt call was made based on a strain that's not in the set, we don't care about that call)

#OG
unique_calls_per_row_OG<- t(apply(OG, 1, function(x) unique(x)))
new_ALT_OG<- lapply(unique_calls_per_row_OG, function(x) paste(x))
new_ALT2_OG<- lapply(new_ALT_OG, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_OG <- data.frame(matrix(unlist(new_ALT2_OG), nrow=length(new_ALT2_OG), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_OG)<- "NEW_ALT_merged"
new_ALT2_OG_sep<- lapply(new_ALT_OG, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_OG_sep <- data.frame(matrix(unlist(new_ALT2_OG_sep), nrow=length(new_ALT2_OG_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_OG_sep)<- "NEW_ALT_OG"

#Clade1
unique_calls_per_row_1<- t(apply(Clade_1, 1, function(x) unique(x)))
new_ALT_1<- lapply(unique_calls_per_row_1, function(x) paste(x))
new_ALT2_1<- lapply(new_ALT_1, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_1 <- data.frame(matrix(unlist(new_ALT2_1), nrow=length(new_ALT2_1), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_1)<- "NEW_ALT_merged"
new_ALT2_1_sep<- lapply(new_ALT_1, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_1_sep <- data.frame(matrix(unlist(new_ALT2_1_sep), nrow=length(new_ALT2_1_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_1_sep)<- "NEW_ALT_Clade_1"

#Clade 2
unique_calls_per_row_2<- t(apply(Clade_2, 1, function(x) unique(x)))
new_ALT_2<- lapply(unique_calls_per_row_2, function(x) paste(x))
new_ALT2_2<- lapply(new_ALT_2, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_2 <- data.frame(matrix(unlist(new_ALT2_2), nrow=length(new_ALT2_2), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_2)<- "NEW_ALT_merged"
new_ALT2_2_sep<- lapply(new_ALT_2, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_2_sep <- data.frame(matrix(unlist(new_ALT2_2_sep), nrow=length(new_ALT2_2_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_2_sep)<- "NEW_ALT_Clade_2"

#all isolates minus clade isolates
unique_calls_per_row_all<- t(apply(All_but_clades, 1, function(x) unique(x)))
new_ALT_all<- lapply(unique_calls_per_row_all, function(x) paste(x))
new_ALT2_all<- lapply(new_ALT_all, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_all <- data.frame(matrix(unlist(new_ALT2_all), nrow=length(new_ALT2_all), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_all)<- "NEW_ALT_merged"
new_ALT2_all_sep<- lapply(new_ALT_all, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_all_sep <- data.frame(matrix(unlist(new_ALT2_all_sep), nrow=length(new_ALT2_all_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_all_sep)<- "NEW_ALT_all"


#add back in the info cols
#OG
new_df_OG<- cbind(OG, new_ALT2_df_OG)
new_df_w_ref_OG<- cbind("NEW_ALT_OG" = new_ALT2_df_OG_sep, new_df_OG, "REF" = snpEff_no_intergenic$REF, "TYPE" = snpEff_no_intergenic$TYPE, "GENE" = snpEff_no_intergenic$GENE, "CHANGEDNA" = snpEff_no_intergenic$CHANGEDNA, "CHANGPEP" = snpEff_no_intergenic$CHANGEPEP, "IMPACT" = snpEff_no_intergenic$IMPACT, POS = snpEff_no_intergenic$POS, CHROM = snpEff_no_intergenic$CHROM)
#C1
new_df_1<- cbind(Clade_1, new_ALT2_df_1)
new_df_w_ref_1<- cbind("NEW_ALT_Clade_1" = new_ALT2_df_1_sep, new_df_1, "REF" = snpEff_no_intergenic$REF, "TYPE" = snpEff_no_intergenic$TYPE, "GENE" = snpEff_no_intergenic$GENE, "CHANGEDNA" = snpEff_no_intergenic$CHANGEDNA, "CHANGPEP" = snpEff_no_intergenic$CHANGEPEP, "IMPACT" = snpEff_no_intergenic$IMPACT, POS = snpEff_no_intergenic$POS, CHROM = snpEff_no_intergenic$CHROM)
#C2
new_df_2<- cbind(Clade_2, new_ALT2_df_2)
new_df_w_ref_2<- cbind("NEW_ALT_Clade_2" = new_ALT2_df_2_sep, new_df_2, "REF" = snpEff_no_intergenic$REF, "TYPE" = snpEff_no_intergenic$TYPE, "GENE" = snpEff_no_intergenic$GENE, "CHANGEDNA" = snpEff_no_intergenic$CHANGEDNA, "CHANGPEP" = snpEff_no_intergenic$CHANGEPEP, "IMPACT" = snpEff_no_intergenic$IMPACT, POS = snpEff_no_intergenic$POS, CHROM = snpEff_no_intergenic$CHROM)
#all
new_df_all<- cbind(All_but_clades, new_ALT2_df_all)
new_df_w_ref_all<- cbind("NEW_ALT_all" = new_ALT2_df_all_sep, new_df_all, "REF" = snpEff_no_intergenic$REF, "TYPE" = snpEff_no_intergenic$TYPE, "GENE" = snpEff_no_intergenic$GENE, "CHANGEDNA" = snpEff_no_intergenic$CHANGEDNA, "CHANGPEP" = snpEff_no_intergenic$CHANGEPEP, "IMPACT" = snpEff_no_intergenic$IMPACT, POS = snpEff_no_intergenic$POS, CHROM = snpEff_no_intergenic$CHROM)

#attach the re-calculated alternative allele collumn from the larger dataset (made by excluding the clades), to the clade datasets
merged_df_C1a<- cbind(NEW_ALT_all = new_df_w_ref_all$NEW_ALT_all, new_df_w_ref_1)
merged_df_C2a<- cbind(NEW_ALT_all = new_df_w_ref_all$NEW_ALT_all, new_df_w_ref_2)

#do the same with the alternative allele collumn calculated using the outgroups
merged_df_C1<- cbind(NEW_ALT_OG = new_df_w_ref_OG$NEW_ALT_OG, merged_df_C1a)
merged_df_C2<- cbind(NEW_ALT_OG = new_df_w_ref_OG$NEW_ALT_OG, merged_df_C2a)



###PROBLEM 1: For each clade, get the number of SNPs presnet in that clade using 1) no targeting, 2) targeting using outgroups, 3) targeting using all isoaltes not in the clades

#subset to only rows where there's a mutation in the clades that's in none of the other isolates i.e. something in the NEW_ALT_CladeX col that's not in the NEW_ALT_all col. 
#clades against all
subset_unique <- function(a,b) FALSE %in%(unlist(strsplit(a,",")) %in% unlist(strsplit(b,",")))
only_in_C1_to_all<- merged_df_C1[apply(merged_df_C1[,c('NEW_ALT_Clade_1','NEW_ALT_all')], 1, function(y) subset_unique(y['NEW_ALT_Clade_1'],y['NEW_ALT_all'])),]
only_in_C2_to_all<- merged_df_C2[apply(merged_df_C2[,c('NEW_ALT_Clade_2','NEW_ALT_all')], 1, function(y) subset_unique(y['NEW_ALT_Clade_2'],y['NEW_ALT_all'])),]

dim(only_in_C1_to_all)
#250
dim(only_in_C2_to_all)  
#87

#clades against out groups
only_in_C1_to_OG<- merged_df_C1[apply(merged_df_C1[,c('NEW_ALT_Clade_1','NEW_ALT_OG')], 1, function(y) subset_unique(y['NEW_ALT_Clade_1'],y['NEW_ALT_OG'])),]
only_in_C2_to_OG<- merged_df_C2[apply(merged_df_C2[,c('NEW_ALT_Clade_2','NEW_ALT_OG')], 1, function(y) subset_unique(y['NEW_ALT_Clade_2'],y['NEW_ALT_OG'])),]

dim(only_in_C1_to_OG)
#3199
dim(only_in_C2_to_OG)  
#643

#accounting for misscalls and dissagreement between the isolates:
only_in_C1_to_OG_confident<- only_in_C1_to_OG[((only_in_C1_to_OG$DMC2_AF100.3B == only_in_C1_to_OG$DMC2_AF100.4B)  & 
                                                 (only_in_C1_to_OG$DMC2_AF100.3B == only_in_C1_to_OG$DMC2_AF100.5B)  & 
                                                 (only_in_C1_to_OG$DMC2_AF100.3B == only_in_C1_to_OG$DMC2_AF100.6B)  & 
                                                 (only_in_C1_to_OG$DMC2_AF100.3B == only_in_C1_to_OG$AF100.10B)  &
                                                 (only_in_C1_to_OG$DMC2_AF100.3B == only_in_C1_to_OG$AF100.12_7)),] 


dim(only_in_C1_to_OG)
dim(only_in_C1_to_OG_confident)

#Clade2
only_in_C2_to_OG_confident<- only_in_C2_to_OG[((only_in_C2_to_OG$DMC2_AF100.8B == only_in_C2_to_OG$DMC2_AF100.9B) & 
                                                 (only_in_C2_to_OG$DMC2_AF100.8B == only_in_C2_to_OG$AF100.10_5)  & 
                                                 (only_in_C2_to_OG$DMC2_AF100.8B == only_in_C2_to_OG$AF100.11_3)  & 
                                                 (only_in_C2_to_OG$DMC2_AF100.8B == only_in_C2_to_OG$DMC2_AF100.12_9)),] 



dim(only_in_C2_to_OG)
dim(only_in_C2_to_OG_confident)

#Deletions
#Clade 1
only_in_C1_to_OG_confident_DELETIONS<- only_in_C1_to_OG[((only_in_C1_to_OG$DMC2_AF100.3B == ".")  & 
                                                           (only_in_C1_to_OG$DMC2_AF100.4B == ".")  & 
                                                           (only_in_C1_to_OG$DMC2_AF100.5B == ".")  & 
                                                           (only_in_C1_to_OG$DMC2_AF100.6B == ".")  & 
                                                           (only_in_C1_to_OG$AF100.10B == ".")  &
                                                           (only_in_C1_to_OG$AF100.12_7 == ".")),] 


dim(only_in_C1_to_OG)
dim(only_in_C1_to_OG_confident_DELETIONS)

#Clade2
only_in_C2_to_OG_confident_DELETIONS<- only_in_C2_to_OG[((only_in_C2_to_OG$DMC2_AF100.8B == ".") & 
                                                           (only_in_C2_to_OG$DMC2_AF100.9B == ".") & 
                                                           (only_in_C2_to_OG$AF100.10_5 == ".")  & 
                                                           (only_in_C2_to_OG$AF100.11_3 == ".")  & 
                                                           (only_in_C2_to_OG$DMC2_AF100.12_9 == ".")),] 

dim(only_in_C2_to_OG)
dim(only_in_C2_to_OG_confident_DELETIONS)

#spanning deletions
#Clade 1
only_in_C1_to_OG_confident_spanning_DELETIONS<- only_in_C1_to_OG[((only_in_C1_to_OG$DMC2_AF100.3B == "*")  & 
                                                                    (only_in_C1_to_OG$DMC2_AF100.4B == "*")  & 
                                                                    (only_in_C1_to_OG$DMC2_AF100.5B == "*")  & 
                                                                    (only_in_C1_to_OG$DMC2_AF100.6B == "*")  & 
                                                                    (only_in_C1_to_OG$AF100.10B == "*")  &
                                                                    (only_in_C1_to_OG$AF100.12_7 == "*")),] 


dim(only_in_C1_to_OG)
dim(only_in_C1_to_OG_confident_spanning_DELETIONS)

#Clade2
only_in_C2_to_OG_confident_spanning_DELETIONS<- only_in_C2_to_OG[((only_in_C2_to_OG$DMC2_AF100.8B == "*") & 
                                                                    (only_in_C2_to_OG$DMC2_AF100.9B == "*") & 
                                                                    (only_in_C2_to_OG$AF100.10_5 == "*")  & 
                                                                    (only_in_C2_to_OG$AF100.11_3 == "*")  & 
                                                                    (only_in_C2_to_OG$DMC2_AF100.12_9 == "*")),] 

dim(only_in_C2_to_OG)
dim(only_in_C2_to_OG_confident_spanning_DELETIONS)


#write to a file
#print results to share
#write.table(only_in_C1_to_OG_confident, "only_in_C1_and_in_OG.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(only_in_C2_to_OG_confident, "only_in_C2_and_in_OG.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##get overlap in variants (all genes with the same variants in C1 and C2, and not in the outgroups)


#do any of these GENES overlap between the two clades? - i.e. found in both clades, in all ioslates but nowhere else in the tree?
#overlaps_by_gene_to_OG<- merge(only_in_C1_to_OG_confident, only_in_C1_to_OG_confident, by.x = "GENE", by.y = "GENE")
overlaps_by_gene_to_OG_specific<- merge(only_in_C1_to_OG_confident, only_in_C2_to_OG_confident, by=c("GENE", "POS"))
overlaps_by_gene_to_OG_genes_only<- merge(only_in_C1_to_OG_confident, only_in_C2_to_OG_confident, by=("GENE"))

#get list of unique genes that intersect
genes_to_pull<- unique(overlaps_by_gene_to_OG_genes_only$GENE)
length(genes_to_pull)
#there are 18 genes in this set

#subset the data frames to include only genes that intersect
genes_C1<- only_in_C1_to_OG_confident[only_in_C1_to_OG_confident$GENE %in% genes_to_pull,]
genes_C2<- only_in_C2_to_OG_confident[only_in_C2_to_OG_confident$GENE %in% genes_to_pull,]

#call table to get number of mutations per gene 
table_genes_C1<- data.frame(table(genes_C1$GENE))
table_genes_C2<- data.frame(table(genes_C2$GENE))

#make sure they're in the same order
table_genes_C1_ordered<-table_genes_C1[order(table_genes_C1$Var1),]
table_genes_C2_ordered<-table_genes_C2[order(table_genes_C2$Var1),]

#combine
print_this<- cbind(table_genes_C1_ordered, table_genes_C2_ordered$Freq )

colnames(print_this)<- c("Gene name", "n variants C1", "n variants C2")

#write.table(print_this, "overlaps_by_gene_to_OG.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# compare this to n SNPS without any subsetting 
#if the REF col matches the new alt, don't return that line (SNP does not exisist in that clade or in the outgroups):
only_in_C1_no_norm<- merged_df_C1[apply(merged_df_C1[,c('NEW_ALT_Clade_1','REF')], 1, function(y) subset_unique(y['NEW_ALT_Clade_1'],y['REF'])),]
only_in_C2_no_norm<- merged_df_C2[apply(merged_df_C2[,c('NEW_ALT_Clade_2','REF')], 1, function(y) subset_unique(y['NEW_ALT_Clade_2'],y['REF'])),]

dim(only_in_C1_no_norm)
#20371
dim(only_in_C2_no_norm)  
#23218


##the above does not account for whether the mutations are in all isolates in a given clade, or just some of them. 
##Here, we look at just mutations that are conserved across each clade. 
#First subset to only SNPS present in all clade isoaltes. 

#subset to only SNPs that appear in all clade isolates but not in the ancesteral strain (ignoring "." unless they're in the ancesteral strain)
#Clade 1
in_all_C1_and_no_others<- only_in_C1_to_all[((only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.4B | only_in_C1_to_all$DMC2_AF100.4B == ".") & 
                                               (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.5B | only_in_C1_to_all$DMC2_AF100.5B == ".") & 
                                               (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.6B | only_in_C1_to_all$DMC2_AF100.6B == ".") & 
                                               (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$AF100.10B | only_in_C1_to_all$AF100.10B == ".") &
                                               (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$AF100.12_7 | only_in_C1_to_all$AF100.12_7 == ".")),] 


dim(only_in_C1_to_all)
dim(in_all_C1_and_no_others)


#Clade 2
in_all_C2_and_no_others<- only_in_C2_to_all[((only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$DMC2_AF100.9B | only_in_C2_to_all$DMC2_AF100.9B == ".") & 
                                               (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$AF100.10_5 | only_in_C2_to_all$AF100.10_5 == ".") & 
                                               (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$AF100.11_3 | only_in_C2_to_all$AF100.11_3 == ".") & 
                                               (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$DMC2_AF100.12_9 | only_in_C2_to_all$DMC2_AF100.12_9 == ".")),] 



dim(only_in_C2_to_all)
dim(in_all_C2_and_no_others)

#the above allows for some isolates to have misscalls/deletions- here, we exclude those calls and only look at true deletions, SNPs, and INDELS
#Clade1
in_all_C1_and_no_others_confident<- only_in_C1_to_all[((only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.4B)  & 
                                                         (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.5B)  & 
                                                         (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$DMC2_AF100.6B)  & 
                                                         (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$AF100.10B)  &
                                                         (only_in_C1_to_all$DMC2_AF100.3B == only_in_C1_to_all$AF100.12_7)),] 


dim(only_in_C1_to_all)
dim(in_all_C1_and_no_others_confident)

#Clade2
in_all_C2_and_no_others_confident<- only_in_C2_to_all[((only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$DMC2_AF100.9B) & 
                                                         (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$AF100.10_5)  & 
                                                         (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$AF100.11_3)  & 
                                                         (only_in_C2_to_all$DMC2_AF100.8B == only_in_C2_to_all$DMC2_AF100.12_9)),] 



dim(only_in_C2_to_all)
dim(in_all_C2_and_no_others_confident)

#print results to share
#write.table(in_all_C1_and_no_others_confident, "only_in_C1_and_in_all.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(in_all_C2_and_no_others_confident, "only_in_C2_and_in_all_C2.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###look for overlaps in unique clade mutations 
#do any of the POSITIONS overlap between the two clades? - i.e. found in both clades, in all ioslates but nowhere else in the tree?
overlaps_by_pos<- merge(in_all_C1_and_no_others_confident, in_all_C2_and_no_others_confident, by.x = "POS", by.y = "POS")
colnames(overlaps_by_pos)<- c("POS", "Clade_1", "Clade_2")
dim(overlaps_by_pos)
#no - there are no overlaps by position

#do any of these GENES overlap between the two clades? - i.e. found in both clades, in all ioslates but nowhere else in the tree?
overlaps_by_gene<- merge(in_all_C1_and_no_others_confident, in_all_C2_and_no_others_confident, by.x = "GENE", by.y = "GENE")
colnames(overlaps_by_gene)<- c("GENE", "Clade2", "Clade_1")
dim(overlaps_by_gene)

#there are two genes that have mutations that don't appear anywhere else in the tree.
#They are both NRPS genes Afu3g03350 and Afu5g12730

#how unexpected is this? How frequently do mutations appear in these two genes in the rest of the tree?
#in the original dataframe (excluding clades), isolate mutations in the genes of interest gene
Afu3g03350<- data.frame(new_df_w_ref_all[new_df_w_ref_all$GENE == "Afu3g03350", ])
Afu5g12730<- data.frame(new_df_w_ref_all[new_df_w_ref_all$GENE == "Afu5g12730", ])

#remove collumn if the call matches the reference
Afu3g03350_subset<- Afu3g03350[,12:length(Afu3g03350) -1]
Afu3g03350_not_in_ref<- data.frame(colSums(Afu3g03350$REF != Afu3g03350_subset))

#remove annotation rows
remove<- c("NEW_ALT_all", 
           "NEW_ALT_merged",
           "REF",
           "TYPE",
           "GENE",
           "CHANGEDNA",
           "CHANGPEP",
           "IMPACT",
           "POS",
           "CHROM")

Afu3g03350_not_in_ref$isolate<- rownames(Afu3g03350_not_in_ref)
Afu3g03350_not_in_ref2<- Afu3g03350_not_in_ref[!(row.names(Afu3g03350_not_in_ref) %in% remove),]
number_of_isolates_with_Afu3g03350_mutations<- sum(Afu3g03350_not_in_ref2$colSums.Afu3g03350.REF....Afu3g03350_subset. > 0)
number_of_isolates_with_Afu3g03350_mutations
#mutated in 80 isoaltes (not including the clades)
max(Afu3g03350_not_in_ref2$colSums.Afu3g03350.REF....Afu3g03350_subset.) 
#from 1 to 12 times 


#remove collumn if the call matches the reference
Afu5g12730_subset<- Afu5g12730[,12:length(Afu5g12730) -1]
Afu5g12730_not_in_ref<- data.frame(colSums(Afu5g12730$REF != Afu5g12730_subset))

#Afu5g12730_not_in_ref2<- Afu5g12730_not_in_ref[ !(row.names(Afu5g12730_not_in_ref) %in% remove), ]
Afu5g12730_not_in_ref$isolate<- rownames(Afu5g12730_not_in_ref)
Afu5g12730_not_in_ref2<- Afu5g12730_not_in_ref[!(row.names(Afu5g12730_not_in_ref) %in% remove),]
number_of_isolates_with_Afu5g12730_mutations<- sum(Afu5g12730_not_in_ref2$colSums.Afu5g12730.REF....Afu5g12730_subset. > 0)
number_of_isolates_with_Afu5g12730_mutations
#mutated in 98 isoaltes (not including the clades)
max(Afu5g12730_not_in_ref2$colSums.Afu5g12730.REF....Afu5g12730_subset.) 
#from 1 to 47 times 
#conclusion - these two genes are frequently mutated. 

#print gene results to share:
unique_genes_clade1<- table(in_all_C1_and_no_others_confident$GENE)
unique_genes_clade2<- table(in_all_C1_and_no_others_confident$GENE)
unique_genes_clade1_ordered<- data.frame(sort(unique_genes_clade1, decreasing = TRUE))
unique_genes_clade2_ordered<- data.frame(sort(unique_genes_clade2, decreasing = TRUE))
colnames(unique_genes_clade1_ordered)<- c("Gene name", "Mutation Freq.")
colnames(unique_genes_clade2_ordered)<- c("Gene name", "Mutation Freq.")

#write.table(unique_genes_clade1_ordered, "only_in_C1_and_in_all_C1_by_gene.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(unique_genes_clade2_ordered, "only_in_C2_and_in_all_C2_by_gene.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




###PROBLEM 2: identify differences between isolates in clade 1 and the ancestral TP1 is0late just outside Clade 1
#C1 ancestral: DMC2_AF100-1_3

#subset to only SNPs that appear in all clade isolates but not in the ancestral strain (ignoring "." unless they're in the ancestral strain)
#with no targeting
all_but_ancesteral_C1<- Clade_1_w_ancestor[((Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.4B | Clade_1_w_ancestor$DMC2_AF100.4B == ".") & 
                                              (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.5B | Clade_1_w_ancestor$DMC2_AF100.5B == ".") & 
                                              (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.6B | Clade_1_w_ancestor$DMC2_AF100.6B == ".") & 
                                              (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$AF100.10B | Clade_1_w_ancestor$AF100.10B == ".") &
                                              (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$AF100.12_7 | Clade_1_w_ancestor$AF100.12_7 == ".") & 
                                              (Clade_1_w_ancestor$DMC2_AF100.3B != Clade_1_w_ancestor$DMC2_AF100.1_3)), ] 


dim(all_but_ancesteral_C1)

#9976 mutations different between the Clade 1 and the ancestral isoalte 

#Excluding possible misscalls and deletions
all_but_ancesteral_C1_confident<- Clade_1_w_ancestor[((Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.4B)  & 
                                                        (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.5B)  & 
                                                        (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$DMC2_AF100.6B)  & 
                                                        (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$AF100.10B)  &
                                                        (Clade_1_w_ancestor$DMC2_AF100.3B == Clade_1_w_ancestor$AF100.12_7) & 
                                                        (Clade_1_w_ancestor$DMC2_AF100.3B != Clade_1_w_ancestor$DMC2_AF100.1_3)),] 

dim(all_but_ancesteral_C1_confident)

#number of confedent deletions 
all_but_ancesteral_C1_confident_DELETIONS<- all_but_ancesteral_C1_confident[((all_but_ancesteral_C1_confident$DMC2_AF100.3B == ".")  & 
                                                                               (all_but_ancesteral_C1_confident$DMC2_AF100.4B == ".")  & 
                                                                               (all_but_ancesteral_C1_confident$DMC2_AF100.5B == ".")  & 
                                                                               (all_but_ancesteral_C1_confident$DMC2_AF100.6B == ".")  &
                                                                               (all_but_ancesteral_C1_confident$AF100.10B == ".") &
                                                                               (all_but_ancesteral_C1_confident$AF100.12_7 == ".")),] 


dim(all_but_ancesteral_C1_confident_DELETIONS)
#there are 2335 deletions the clade 1 relative to TP1


#Targeting with the outgroups prior to analysis:
#bind all relevant rows from TP1 to only_in_C1_to_OG
C1_to_OG_w_ancestor<- cbind(only_in_C1_to_OG, TP1 = snpEff_no_intergenic[, "DMC2_AF100.1_3"][match(only_in_C1_to_OG$POS, snpEff_no_intergenic$POS)])

#subset to only those rows matching 
all_but_ancesteral_C1_using_OG<- C1_to_OG_w_ancestor[((C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.4B | C1_to_OG_w_ancestor$DMC2_AF100.4B == ".") & 
                                                        (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.5B | C1_to_OG_w_ancestor$DMC2_AF100.5B == ".") & 
                                                        (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.6B | C1_to_OG_w_ancestor$DMC2_AF100.6B == ".") & 
                                                        (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$AF100.10B | C1_to_OG_w_ancestor$AF100.10B == ".") &
                                                        (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$AF100.12_7 | C1_to_OG_w_ancestor$AF100.12_7 == ".") & 
                                                        (C1_to_OG_w_ancestor$DMC2_AF100.3B != C1_to_OG_w_ancestor$TP1)), ] 


dim(C1_to_OG_w_ancestor)
dim(all_but_ancesteral_C1_using_OG)

#Excluding possible misscalls and deletions
all_but_ancesteral_C1_confident_using_OG<- C1_to_OG_w_ancestor[((C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.4B)  & 
                                                                   (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.5B)  & 
                                                                   (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$DMC2_AF100.6B)  & 
                                                                   (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$AF100.10B)  &
                                                                   (C1_to_OG_w_ancestor$DMC2_AF100.3B == C1_to_OG_w_ancestor$AF100.12_7) & 
                                                                   (C1_to_OG_w_ancestor$DMC2_AF100.3B != C1_to_OG_w_ancestor$TP1)),] 

dim(C1_to_OG_w_ancestor)
dim(all_but_ancesteral_C1_confident_using_OG)

#number of confident deletions 
all_but_ancesteral_C1_confident_using_OG_DELETIONS<- all_but_ancesteral_C1_confident_using_OG[((all_but_ancesteral_C1_confident_using_OG$DMC2_AF100.3B == ".")  & 
                                                                                                 (all_but_ancesteral_C1_confident_using_OG$DMC2_AF100.4B == ".")  & 
                                                                                                 (all_but_ancesteral_C1_confident_using_OG$DMC2_AF100.5B == ".")  & 
                                                                                                 (all_but_ancesteral_C1_confident_using_OG$DMC2_AF100.6B == ".")  &
                                                                                                 (all_but_ancesteral_C1_confident_using_OG$AF100.10B == ".") &
                                                                                                 (all_but_ancesteral_C1_confident_using_OG$AF100.12_7 == ".")),] 


dim(all_but_ancesteral_C1_confident_using_OG_DELETIONS)


##look by gene
#call table to get number of mutations in each gene
genes_mutated_in_C1_but_not_in_ancestor<- table(all_but_ancesteral_C1_confident_using_OG$GENE)

genes_mutated_in_C1_but_not_in_ancestor_ordered<- data.frame(sort(genes_mutated_in_C1_but_not_in_ancestor, decreasing = TRUE))
colnames(genes_mutated_in_C1_but_not_in_ancestor_ordered)<- c("Gene name", "Mutation Freq.")


#print to file
#write.table(genes_mutated_in_C1_but_not_in_ancestor_ordered, "genes_mutated_in_C1_but_not_in_ancestor.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(all_but_ancesteral_C1_confident_using_OG, "mutations_in_C1_but_not_in_ancestor.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



###Problem 3:
#Look at mutational acumulation over time. 
# first get synonamous mutations to compare against 
snpEff_synonymous<- no_intergenic_b[no_intergenic_b$TYPE == "synonymous_variant",]

#subset the clades from the main df (for syn mutations)
OG_synonymous<- snpEff_synonymous[, OG_names]
Clade_1_synonymous<- snpEff_synonymous[, Clade_1_names]
Clade_2_synonymous<- snpEff_synonymous[, Clade_2_names]

#target intergenic variants using outgroups
#OG
unique_calls_per_row_OG<- t(apply(OG_synonymous, 1, function(x) unique(x)))
new_ALT_OG<- lapply(unique_calls_per_row_OG, function(x) paste(x))
new_ALT2_OG<- lapply(new_ALT_OG, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_OG <- data.frame(matrix(unlist(new_ALT2_OG), nrow=length(new_ALT2_OG), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_OG)<- "NEW_ALT_merged"
new_ALT2_OG_sep<- lapply(new_ALT_OG, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_OG_sep_syn <- data.frame(matrix(unlist(new_ALT2_OG_sep), nrow=length(new_ALT2_OG_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_OG_sep_syn)<- "NEW_ALT_OG"

#Clade1
unique_calls_per_row_1<- t(apply(Clade_1_synonymous, 1, function(x) unique(x)))
new_ALT_1<- lapply(unique_calls_per_row_1, function(x) paste(x))
new_ALT2_1<- lapply(new_ALT_1, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_1 <- data.frame(matrix(unlist(new_ALT2_1), nrow=length(new_ALT2_1), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_1)<- "NEW_ALT_merged"
new_ALT2_1_sep<- lapply(new_ALT_1, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_1_sep_syn <- data.frame(matrix(unlist(new_ALT2_1_sep), nrow=length(new_ALT2_1_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_1_sep_syn)<- "NEW_ALT_Clade_1"

#Clade 2
unique_calls_per_row_2<- t(apply(Clade_2_synonymous, 1, function(x) unique(x)))
new_ALT_2<- lapply(unique_calls_per_row_2, function(x) paste(x))
new_ALT2_2<- lapply(new_ALT_2, function(x) paste(x, sep = " ", collapse = ""))
new_ALT2_df_2 <- data.frame(matrix(unlist(new_ALT2_2), nrow=length(new_ALT2_2), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_2)<- "NEW_ALT_merged"
new_ALT2_2_sep<- lapply(new_ALT_2, function(x) paste(x, sep = " ", collapse = ","))
new_ALT2_df_2_sep_syn <- data.frame(matrix(unlist(new_ALT2_2_sep), nrow=length(new_ALT2_2_sep), byrow=T),stringsAsFactors=FALSE)
colnames(new_ALT2_df_2_sep_syn)<- "NEW_ALT_Clade_2"

#add back in the info cols
new_df_OG<- cbind(OG_synonymous, new_ALT2_df_OG)
new_df_w_ref_OG<- cbind("NEW_ALT_OG" = new_ALT2_df_OG_sep_syn, new_df_OG, "REF" = snpEff_synonymous$REF, "TYPE" = snpEff_synonymous$TYPE, "GENE" = snpEff_synonymous$GENE, "CHANGEDNA" = snpEff_synonymous$CHANGEDNA, "CHANGPEP" = snpEff_synonymous$CHANGEPEP, "IMPACT" = snpEff_synonymous$IMPACT, POS = snpEff_synonymous$POS, CHROM = snpEff_synonymous$CHROM)
#C1
new_df_1<- cbind(Clade_1_synonymous, new_ALT2_df_1)
new_df_w_ref_1<- cbind("NEW_ALT_Clade_1" = new_ALT2_df_1_sep_syn, new_df_1, "REF" = snpEff_synonymous$REF, "TYPE" = snpEff_synonymous$TYPE, "GENE" = snpEff_synonymous$GENE, "CHANGEDNA" = snpEff_synonymous$CHANGEDNA, "CHANGPEP" = snpEff_synonymous$CHANGEPEP, "IMPACT" = snpEff_synonymous$IMPACT, POS = snpEff_synonymous$POS, CHROM = snpEff_synonymous$CHROM)
#C2
new_df_2<- cbind(Clade_2_synonymous, new_ALT2_df_2)
new_df_w_ref_2<- cbind("NEW_ALT_Clade_2" = new_ALT2_df_2_sep_syn, new_df_2, "REF" = snpEff_synonymous$REF, "TYPE" = snpEff_synonymous$TYPE, "GENE" = snpEff_synonymous$GENE, "CHANGEDNA" = snpEff_synonymous$CHANGEDNA, "CHANGPEP" = snpEff_synonymous$CHANGEPEP, "IMPACT" = snpEff_synonymous$IMPACT, POS = snpEff_synonymous$POS, CHROM = snpEff_synonymous$CHROM)

#do the same with the alternative allele collumn calculated using the outgroups
merged_df_C1_syn<- cbind(NEW_ALT_OG = new_df_w_ref_OG$NEW_ALT_OG, new_df_w_ref_1)
merged_df_C2_syn<- cbind(NEW_ALT_OG = new_df_w_ref_OG$NEW_ALT_OG, new_df_w_ref_2)

dim(merged_df_C1_syn)
dim(merged_df_C2_syn)

#subset based on OG targeting
#clades against out groups
only_in_C1_to_OG_syn<- merged_df_C1_syn[apply(merged_df_C1_syn[,c('NEW_ALT_Clade_1','NEW_ALT_OG')], 1, function(y) subset_unique(y['NEW_ALT_Clade_1'],y['NEW_ALT_OG'])),]
only_in_C2_to_OG_syn<- merged_df_C2_syn[apply(merged_df_C2_syn[,c('NEW_ALT_Clade_2','NEW_ALT_OG')], 1, function(y) subset_unique(y['NEW_ALT_Clade_2'],y['NEW_ALT_OG'])),]

dim(only_in_C1_to_OG_syn)
dim(only_in_C2_to_OG_syn)

dim(only_in_C1_to_OG)
dim(only_in_C2_to_OG)


#is this reflective of starting number?  
dim(snpEff_synonymous)
dim(snpEff_no_intergenic)
#yes

#ACUMULATION ANALYSIS
####Step 1
#subsets the input dfs to only include the strains in each clade
convert_to_t_f <- function(df, timepts){
  subsetdf<- data.frame(df[, grepl(paste(timepts, collapse="|"),names(df)) ])
  T_F_df <- data.frame(df$REF == subsetdf[,1: ncol(subsetdf)])
  return(T_F_df)
}


#set timepoints
timepts_clade1<- as.numeric(c(3,4,5,6,10,12)) 
timepts_clade2<- as.numeric(c(8,9,10,11,12))

#run functon
T_F_df_clade1<- convert_to_t_f(df = only_in_C1_to_OG, timepts = timepts_clade1)
T_F_df_clade2<- convert_to_t_f(df = only_in_C2_to_OG, timepts = timepts_clade2)

T_F_df_clade1_syn<- convert_to_t_f(df = only_in_C1_to_OG_syn, timepts = timepts_clade1)
T_F_df_clade2_syn<- convert_to_t_f(df = only_in_C2_to_OG_syn, timepts = timepts_clade2)


####Step 2
#look, only at the SNPS that have gone to fixation
#define function 
Rcpp::cppFunction('
CharacterVector all_T(LogicalMatrix x) {
  CharacterVector col_names = colnames(x);
  CharacterVector out(x.nrow(), col_names(0));
  for (int i = 0; i < x.nrow(); i++){
    for (int j = x.ncol() - 1; j >= 0; j--){
      if (!x(i, j)){
        if (j == x.ncol() - 1) {
          out(i) = NA_STRING;
          } else {
            out(i) = col_names(j + 1);
          }
          break;
          }
      }
    }
  return(out);
}
                  ')

fixed_muts_appeared_at_what_TP_C1<- data.frame(all_T(as.matrix(T_F_df_clade1)))
fixed_muts_appeared_at_what_TP_C2<- data.frame(all_T(as.matrix(T_F_df_clade2)))
row.names(fixed_muts_appeared_at_what_TP_C1)<- rownames(T_F_df_clade1)
row.names(fixed_muts_appeared_at_what_TP_C2)<- rownames(T_F_df_clade2)

fixed_muts_appeared_at_what_TP_C1_syn<- data.frame(all_T(as.matrix(T_F_df_clade1_syn)))
fixed_muts_appeared_at_what_TP_C2_syn<- data.frame(all_T(as.matrix(T_F_df_clade2_syn)))
row.names(fixed_muts_appeared_at_what_TP_C1_syn)<- rownames(T_F_df_clade1_syn)
row.names(fixed_muts_appeared_at_what_TP_C2_syn)<- rownames(T_F_df_clade2_syn)


#call table to get totals
table_C1<- data.frame(table(fixed_muts_appeared_at_what_TP_C1))
table_C2<- data.frame(table(fixed_muts_appeared_at_what_TP_C2))

table_C1_syn<- data.frame(table(fixed_muts_appeared_at_what_TP_C1_syn))
table_C2_syn<- data.frame(table(fixed_muts_appeared_at_what_TP_C2_syn))


#order the tables 
order_C1 <- c("DMC2_AF100.3B", 
              "DMC2_AF100.4B", 
              "DMC2_AF100.5B", 
              "DMC2_AF100.6B",
              "AF100.10B", 
              "AF100.12_7")

order_C2 <- c("DMC2_AF100.8B", 
              "DMC2_AF100.9B", 
              "AF100.10_5", 
              "AF100.11_3",
              "DMC2_AF100.12_9")

table_C1_syn
#C1_syn is zero for TPs 4 and 6
table_C2_syn
#C2_syn is zero for TP 8
#add them back in 
add_back_C1<- data.frame(rbind(cbind("DMC2_AF100.6B", 0), cbind("DMC2_AF100.4B", 0)))
colnames(add_back_C1)<- c("fixed_muts_appeared_at_what_TP_C1_syn", "Freq")
add_back_C2<- data.frame(cbind("DMC2_AF100.8B", 0))
colnames(add_back_C2)<- c("fixed_muts_appeared_at_what_TP_C2_syn", "Freq")
table_C1_syn<- rbind(table_C1_syn, add_back_C1)
table_C2_syn<- rbind(table_C2_syn, add_back_C2)


table_C1_ordered <- table_C1 %>% slice(match(order_C1,fixed_muts_appeared_at_what_TP_C1))
table_C2_ordered <- table_C2 %>% slice(match(order_C2,fixed_muts_appeared_at_what_TP_C2))

table_C1_ordered_syn <- table_C1_syn %>% slice(match(order_C1,fixed_muts_appeared_at_what_TP_C1_syn))
table_C2_ordered_syn <- table_C2_syn %>% slice(match(order_C2,fixed_muts_appeared_at_what_TP_C2_syn))

#check
table_C1_ordered 
table_C2_ordered

#attach time pts (days since ioslation)
dsi1<- c(135, 249, 471, 624, 1221, 1498)
dsi2<- c(968, 1145, 1221, 1346, 1498)

table_C1_ordered$dsi<- dsi1
table_C2_ordered$dsi<- dsi2

#attach synonamous mutaitons 
table_C1_ordered$syn_muts<- table_C1_ordered_syn$Freq
table_C2_ordered$syn_muts<- table_C2_ordered_syn$Freq


##GRAPH fixed mutations over time (excluding first T pt.)
#Clade1
plot(table_C1_ordered$dsi[2:nrow(table_C1_ordered)], table_C1_ordered$Freq[2:nrow(table_C1_ordered)],
     col = "#50184499", 
     xlab = "days after isolation", 
     ylab = "total new mutations",
     main = "Clade 1 fixed mutations", 
     pch = 16, 
     xlim=c(0,1500),
     ylim=c(0,120))

clade1_new_lm<- lm(table_C1_ordered$Freq[2:nrow(table_C1_ordered)]~ table_C1_ordered$dsi[2:nrow(table_C1_ordered)])
clade1_new_lm_zero<- lm(table_C1_ordered$Freq[2:nrow(table_C1_ordered)]~ 0 + table_C1_ordered$dsi[2:nrow(table_C1_ordered)])
sum1<- summary(clade1_new_lm)
sum1_C1_forced_to_zero<- summary(clade1_new_lm_zero)
anova(clade1_new_lm, clade1_new_lm_zero)
#not significantly different, but RSS is higher for the model with zero thoguh the origin, so use that one. 
#slope for that model is 0.05355
1/0.05355

abline(clade1_new_lm, col = "#501844")
rsq_clade1_new <- bquote(italic(R)^2 == .(format(summary(clade1_new_lm)$adj.r.squared, digits=2)))
text(x = 200, y = 95, labels = rsq_clade1_new, col = "#501844", cex = .8)
pval_calde1_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde1_new <- bquote(italic(P) == .(pval_calde1_new))
text(x = 200, y = 110, labels = pval_calde1_new, col = "#501844", cex = .8)
#mutations are increasing over time, but not significantly in Calde 1
par(new=T)

plot(table_C1_ordered$dsi[2:nrow(table_C1_ordered)], table_C1_ordered$syn_muts[2:nrow(table_C1_ordered)],
     col = "#F0593999", 
     xlab = "", 
     ylab = "",
     main = "", 
     axes=F,
     pch = 16,
     xlim=c(0,1500),
     ylim=c(0,120))
clade1_new_lm<- lm(table_C1_ordered$syn_muts[2:nrow(table_C1_ordered)]~ table_C1_ordered$dsi[2:nrow(table_C1_ordered)])
sum1<- summary(clade1_new_lm)
abline(clade1_new_lm, col = "#F05939")
rsq_clade1_new <- bquote(italic(R)^2 == .(format(summary(clade1_new_lm)$adj.r.squared, digits=2)))
text(x = 200, y = 55, labels = rsq_clade1_new, col = "#F05939", cex = .8)
pval_calde1_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde1_new <- bquote(italic(P) == .(pval_calde1_new))
text(x = 200, y = 70, labels = pval_calde1_new, col = "#F05939", cex = .8)


par(opar)  # restore plot par

#Clade 2
plot(table_C2_ordered$dsi[2:nrow(table_C2_ordered)], table_C2_ordered$Freq[2:nrow(table_C2_ordered)],
     col = "#50184499", 
     xlab = "days after isolation", 
     ylab = "total new mutations",
     main = "Clade 2 fixed mutations", 
     pch = 16, 
     xlim=c(1100,1500),
     ylim=c(0,120))

clade2_new_lm<- lm(table_C2_ordered$Freq[2:nrow(table_C2_ordered)]~ table_C2_ordered$dsi[2:nrow(table_C2_ordered)])
#to report slope - find out if you need to force the intercept to 0 = run one version 
#forced to zero and see if one has a significantly higher RSS value - use that one. 
clade2_new_lm_zero<- lm(table_C2_ordered$Freq[2:nrow(table_C2_ordered)]~ 0 + table_C2_ordered$dsi[2:nrow(table_C2_ordered)])
sum1<- summary(clade2_new_lm)
sum1_C2_forced_to_zero<- summary(clade2_new_lm_zero)
anova(clade2_new_lm, clade2_new_lm_zero)
#model with intercept forced to zero is better. slope ("Estimate") is 0.03799
#to get the number of days per one mutation fixed:
1/0.03799

abline(clade2_new_lm, col = "#501844")
rsq_clade2_new <- bquote(italic(R)^2 == .(format(summary(clade2_new_lm)$adj.r.squared, digits=2)))
text(x = 1150, y = 95, labels = rsq_clade2_new, col = "#501844", cex = .8)
pval_calde2_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde2_new <- bquote(italic(P) == .(pval_calde2_new))
text(x = 1150, y = 110, labels = pval_calde2_new, col = "#501844", cex = .8)
#mutations are increasing over time, but not significantly in Calde 1
par(new=T)

plot(table_C2_ordered$dsi[2:nrow(table_C2_ordered)], table_C2_ordered$syn_muts[2:nrow(table_C2_ordered)],
     col = "#F0593999", 
     xlab = "", 
     ylab = "",
     main = "", 
     axes=F,
     pch = 16,
     xlim=c(1100,1500),
     ylim=c(0,120))
clade1_new_lm<- lm(table_C2_ordered$syn_muts[2:nrow(table_C2_ordered)]~ table_C2_ordered$dsi[2:nrow(table_C2_ordered)])
sum1<- summary(clade1_new_lm)
abline(clade1_new_lm, col = "#F05939")
rsq_clade1_new <- bquote(italic(R)^2 == .(format(summary(clade1_new_lm)$adj.r.squared, digits=2)))
text(x = 1150, y = 55, labels = rsq_clade1_new, col = "#F05939", cex = .8)
pval_calde1_new<- round(sum1$coefficients[2,4], digits = 3)
pval_calde1_new <- bquote(italic(P) == .(pval_calde1_new))
text(x = 1150, y = 70, labels = pval_calde1_new, col = "#F05939", cex = .8)
#non-synonymous mutations are increasing significantly over time in Calde 2
#add legend
legend("topleft", inset=.02,c("synonymous","non-synonymous"), fill=c("#F0593999", "#50184499"), horiz=TRUE, cex=0.8)



##subset files to print of the mutations that have gone to fixation along with the time point at which they became fixed
#use rownames from fixed_muts_appeared_at_what_TP_C1 and attach the relevent collumns from 
all_info_C1<- cbind(only_in_C1_to_OG, When_fixed_mut_appeared = fixed_muts_appeared_at_what_TP_C1$all_T.as.matrix.T_F_df_clade1..)
all_info_C2<- cbind(only_in_C2_to_OG, When_fixed_mut_appeared = fixed_muts_appeared_at_what_TP_C2$all_T.as.matrix.T_F_df_clade2..)

#subset to only fixed SNPs 
all_info_C1_fixed<- all_info_C1[!is.na(all_info_C1$When_fixed_mut_appeared),]
all_info_C2_fixed<- all_info_C2[!is.na(all_info_C2$When_fixed_mut_appeared),]

#check that dimensions are correct
dim(all_info_C1_fixed)
dim(all_info_C2_fixed)
#compared to:
sum(table_C1_ordered$Freq)
sum(table_C2_ordered$Freq)
#for synonamous 
sum(as.numeric(table_C1_ordered$syn_muts))
sum(as.numeric(table_C2_ordered$syn_muts))

#print these to share:
#write.table(all_info_C1_fixed, "fixed_mutations_C1.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#write.table(all_info_C2_fixed, "fixed_mutations_C2.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


###CALCULATE RATES
##get mutation rate (need to calculate the number of days since last time point)
table_C1_ordered$days1<- c(135, 249 - 135, 471 - 249, 624 - 471, 1221 - 624, 1498 - 1221)
table_C1_ordered$rate<- table_C1_ordered$Freq / table_C1_ordered$days1
fixed_totals_per_timept1 <- table_C1_ordered[2:6,]
muts_fixed_per_day_C1<- mean(fixed_totals_per_timept1$rate)
non_syn_fixed_C1<- 1/(muts_fixed_per_day_C1) 


#for clade 2
table_C2_ordered$days1<- c(968, 1145 - 968, 1221 - 1145, 1346 - 1221, 1498 - 1346)
table_C2_ordered$rate<- table_C2_ordered$Freq / table_C2_ordered$days1
fixed_totals_per_timept2 <- table_C2_ordered[2:5,]
muts_fixed_per_day_C2<- mean(fixed_totals_per_timept2$rate)
non_syn_fixed_C2<- 1/(muts_fixed_per_day_C2) 

