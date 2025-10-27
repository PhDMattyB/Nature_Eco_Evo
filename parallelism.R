######################################
## genomic parallelism test
##
## 08.10.2025
##
## Matthew K Brachmann MKB phdmattyb
######################################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

library(tidyverse)
library(ggVennDiagram)

## (Sharedloci ∕ totallociinpopulation 1 + Sharedloci ∕ totallociinpopulation 2) ∗ 0.5

# Fst outlier loci --------------------------------------------------------

ASHN_out_loc = read_csv('ASHN_TOP_DAWG_FST_outlier.csv') %>% 
  dplyr::select(SNP)
MYV_out_loc = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP)
SKR_out_loc = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP)
GTS_GAR_out_loc = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP)

list_out_loc = list(ASHN_out_loc, 
     MYV_out_loc, 
     SKR_out_loc, 
     GTS_GAR_out_loc)

ggVennDiagram(list_out_loc)

# parallelism per loci ----------------------------------------------------

## ASHN vs myv
((0/281)+(0/267)) * 0.5

## ASH vs skr
((1/281)+(1/292)) * 0.5

## ASHN vs gts-gar
((1/281)+(1/351)) * 0.5

## MYV vs skr
((1/267)+ (1/292))*0.5

## myv vs gts-gar
((1/267)+ (1/351))*0.5

## SKR vs gts-gar
((3/292) + (3/351))*0.5



# gene parallelims --------------------------------------------------------


ASHN_genes = read_tsv('ASHN_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                      col_names = F) %>% 
  distinct(X1) %>% 
  separate(X1, 
           c('gene_name', 
             'trash'),
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)
MYV_genes = read_tsv("MYV_FST_0.5%_outlier_gene_names_only_NEW.tsv", 
                     col_names = F)%>% 
  distinct(X1) %>% 
  separate(X1, 
           c('gene_name', 
             'trash'),
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)
SKR_genes = read_tsv("SKR_FST_0.5%_outlier_gene_names_only_NEW.tsv", 
                     col_names = F)%>% 
  distinct(X1) %>% 
  separate(X1, 
           c('gene_name', 
             'trash'),
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)
GTS_gar_genes = read_tsv('GTS_CSWY_FST_0.5%_outlier_gene_names_only_NEW.tsv', 
                         col_names = F)%>% 
  distinct(X1) %>% 
  separate(X1, 
           c('gene_name', 
             'trash'),
           sep = '-') %>% 
  dplyr::select(-trash) %>% 
  distinct(gene_name)

## ASHN vs myv
((0/401)+(0/393)) * 0.5

## ASH vs skr
((1/401)+(1/434)) * 0.5

## ASHN vs gts-gar
((1/401)+(1/460)) * 0.5

## MYV vs skr
((4/393)+ (4/434))*0.5

## myv vs gts-gar
((2/393)+ (2/460))*0.5

## SKR vs gts-gar
((3/434) + (3/460))*0.5




# Sorensen dice similarity index - loci ------------------------------------------

sorensen_dice <- (2 * shared_loci) / (total_pop1 + total_pop2)


## ASHN vs myv
(2*0)/(281+267)

## ASH vs skr
(2*1)/(281+292)

## ASHN vs gts-gar
(2*1)/(281+351)

## MYV vs skr
(2*1)/(267+292)
## myv vs gts-gar
(2*2)/(267+351)
## SKR vs gts-gar
(2*3)/(282+351)
# fisher exact test - loci ------------------------------------------------

ASHN_MYV_loci = matrix(c(0, 281, 267, 0), 
                       nrow = 2, 
                       byrow = T)

fisher.test(ASHN_MYV_loci)

ASHN_skr_loci = matrix(c(1, 281, 292, 0), 
                       nrow = 2, 
                       byrow = T)

fisher.test(ASHN_skr_loci)

ASHN_gtsgar_loci = matrix(c(1, 281, 351, 0), 
                       nrow = 2, 
                       byrow = T)

fisher.test(ASHN_gtsgar_loci)

myv_skr_loci = matrix(c(1, 267, 292, 0), 
                         nrow = 2, 
                         byrow = T)

fisher.test(myv_skr_loci)

myv_gtsgar_loci = matrix(c(2, 267, 351, 0), 
                          nrow = 2, 
                          byrow = T)

fisher.test(myv_gtsgar_loci)

skr_gtsgar_loci = matrix(c(3, 292, 351, 0), 
                      nrow = 2, 
                      byrow = T)

fisher.test(skr_gtsgar_loci)



# Sorensen dice similarity index - Genes ------------------------------------------

##ASHN vs myv
(2*0)/(401+393)

## ASH vs skr
(2*1)/(401+434)
## ASHN vs gts-gar
(2*1)/(401+460)

## MYV vs skr
(2*4)/(393+434)
## myv vs gts-gar
(2*2)/(393+460)
## SKR vs gts-gar
(2*3)/(434+460)


# fisher exact genes  -----------------------------------------------------

ASHN_MYV_gene = matrix(c(0, 401, 393, 0), 
                       nrow = 2, 
                       byrow = T)

fisher.test(ASHN_MYV_gene)

ASHN_skr_gene = matrix(c(1, 401, 434, 0), 
                       nrow = 2, 
                       byrow = T)

fisher.test(ASHN_skr_gene)

ASHN_gtsgar_gene = matrix(c(1, 401, 460, 0), 
                          nrow = 2, 
                          byrow = T)

fisher.test(ASHN_gtsgar_gene)

myv_skr_gene = matrix(c(4, 393, 434, 0), 
                      nrow = 2, 
                      byrow = T)

fisher.test(myv_skr_gene)

myv_gtsgar_gene = matrix(c(2, 393, 460, 0), 
                         nrow = 2, 
                         byrow = T)

fisher.test(myv_gtsgar_gene)

skr_gtsgar_gene = matrix(c(3, 434, 460, 0), 
                         nrow = 2, 
                         byrow = T)

fisher.test(skr_gtsgar_gene)

