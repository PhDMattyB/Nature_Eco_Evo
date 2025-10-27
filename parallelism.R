######################################
## genomic parallelism test
##
## 08.10.2025
##
## Matthew K Brachmann MKB phdmattyb
######################################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

library(tidyverse)
library(ggvenn)
library(ggVennDiagram)

## (Sharedloci ∕ totallociinpopulation 1 + Sharedloci ∕ totallociinpopulation 2) ∗ 0.5

# Fst outlier loci --------------------------------------------------------

ASHN_out_loc = read_csv('ASHN_TOP_DAWG_FST_outlier.csv') %>% 
  dplyr::select(SNP) %>% 
  as.list()
MYV_out_loc = read_csv('MYV_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP) %>% 
  as.list()
SKR_out_loc = read_csv('SKR_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP) %>% 
  as.list()
GTS_GAR_out_loc = read_csv('GTS_CSWY_TOP_DAWG_Fst_clean.csv') %>% 
  filter(value == 'Outlier') %>% 
  dplyr::select(SNP) %>% 
  as.list()

list_out_loc = list(ASHN_out_loc$SNP, 
     MYV_out_loc$SNP, 
     SKR_out_loc$SNP, 
     GTS_GAR_out_loc$SNP)

location_cols = c('#00798c',
                  '#003d5b',
                  '#edae49',
                  '#d1495b',
                  '#30638e')

out_loc_venn = ggVennDiagram(list_out_loc, 
              label_alpha = 0,
              label = 'count',
              category.names = c("ASHN - Young",
                                 "MYV - Old",
                                 "SKR - Young", 
                                 'GTS/GAR - Old'),
              label_size = 4, 
              set_color = c("#00798c",
                            "#d1495b",
                            "#30638e",
                            "#edae49"))+
  scale_fill_gradient(low = "#ffffff", high = "#ffffff")+ 
  scale_x_continuous(expand = expansion(mult = .2))+
  theme_void()+
  theme(
  text = element_text(size = 14), 
  legend.position = 'none')+
  labs(title = "F) Parallelism of outlier loci")


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


ASHN_genes = read_csv('ASHN_FST_0.5%_outlier_genes_NEW.csv', 
                      col_names = T) %>% 
  filter(feature == 'gene') %>% 
  filter(., !grepl('ENSG', gene_name)) %>% 
  arrange(gene_name) %>% 
  distinct(gene_name)
# %>% 
#   distinct(X1) %>% 
#   separate(X1, 
#            c('gene_name', 
#              'trash'),
#            sep = '-') %>% 
#   dplyr::select(-trash) %>% 
#   distinct(gene_name)
MYV_genes = read_csv("MYV_FST_0.5%_outlier_genes_NEW.csv", 
                     col_names = T)%>% 
  filter(feature == 'gene') %>% 
  filter(., !grepl('ENSG', gene_name)) %>% 
  arrange(gene_name) %>% 
  distinct(gene_name)
  # distinct(X1) %>% 
  # separate(X1, 
  #          c('gene_name', 
  #            'trash'),
  #          sep = '-') %>% 
  # dplyr::select(-trash) %>% 
  # distinct(gene_name)
SKR_genes = read_csv("SKR_FST_0.5%_outlier_genes_NEW.csv", 
                     col_names = T)%>% 
  filter(feature == 'gene') %>% 
  filter(., !grepl('ENSG', gene_name)) %>% 
  arrange(gene_name) %>% 
  distinct(gene_name)
  
GTS_gar_genes = read_csv('GTS_CSWY_FST_0.5%_outlier_genes_NEW.csv', 
                         col_names = T)%>% 
  filter(feature == 'gene') %>% 
  filter(., !grepl('ENSG', gene_name)) %>% 
  arrange(gene_name) %>% 
  distinct(gene_name)


list_out_gene = list(ASHN_genes$gene_name, 
                    MYV_genes$gene_name, 
                    SKR_genes$gene_name, 
                    GTS_gar_genes$gene_name)


out_gene_venn = ggVennDiagram(list_out_gene, 
                             label_alpha = 0,
                             label = 'count',
                             category.names = c("ASHN - Young",
                                                "MYV - Old",
                                                "SKR - Young", 
                                                'GTS/GAR - Old'),
                             label_size = 4, 
                             set_color = c("#00798c",
                                           "#d1495b",
                                           "#30638e",
                                           "#edae49"))+
  scale_fill_gradient(low = "#ffffff", high = "#ffffff")+ 
  scale_x_continuous(expand = expansion(mult = .2))+
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = 'none')+
  labs(title = "G) Parallelism of genes")





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


# simulated example data --------------------------------------------------

sim_data = read_csv("VENN_simulated_data.csv")

sim_ASHN = sim_data %>% 
  filter(ecotype_pair == 'ASHN') %>% 
  dplyr::select(SNP)
sim_SKR = sim_data %>% 
  filter(ecotype_pair == 'SKR')%>% 
  dplyr::select(SNP)
sim_MYV = sim_data %>% 
  filter(ecotype_pair == 'MYV')%>% 
  dplyr::select(SNP)
sim_GTS = sim_data %>% 
  filter(ecotype_pair == 'GTS')%>% 
  dplyr::select(SNP)


# intersect(sim_GTS,
#           sim_MYV) %>%
#   intersect(.,
#             sim_GTS) %>%
#   intersect(.,
#             sim_GTS)

list_out_sim = list(sim_ASHN$SNP, 
                     sim_MYV$SNP, 
                     sim_SKR$SNP, 
                     sim_GTS$SNP)

out_sim_venn = ggVennDiagram(list_out_sim, 
              label_alpha = 0,
              label = 'count',
              category.names = c("ASHN - Young",
                                 "MYV - Old",
                                 "SKR - Young", 
                                 'GTS/GAR - Old'),
              label_size = 4, 
              set_color = c("#00798c",
                            "#d1495b",
                            "#30638e",
                            "#edae49"))+
  scale_fill_gradient(low = "#ffffff", high = "#ffffff")+ 
  scale_x_continuous(expand = expansion(mult = .2))+
  theme_void()+
  theme(
    text = element_text(size = 14), 
    legend.position = 'none')+
  labs(title = "E) Predicted level of parallelism")


# combine venn diagrams ---------------------------------------------------

venn_combo = out_sim_venn/out_loc_venn/out_gene_venn
venn_combo = out_sim_venn + out_loc_venn + out_gene_venn

ggsave('VennDiagram_Genomic_Parallelism.svg', 
       plot = venn_combo, 
       dpi = 'retina', 
       units = 'cm', 
       width = 45, 
       height = 15)

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

