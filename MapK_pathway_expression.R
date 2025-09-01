##############################
## MAPK signalling divergence
##
## Matt Brachmann (PhDMattyB)
##
## 01.09.2025
##
##############################


library(tidyverse)

setwd('~/Parsons_Postdoc/Stickleback_Genomic/Stickleback_Annotation_features/')

brain_f1_exp = read_csv('BRAIN_BH_FDR_GLMER_ecotemp_pval0.05.csv')


liver_f1_exp = read_csv('LIVER_BH_FDR_GLMER_ecotemp_pval0.05.csv')

# brain_f1_exp %>% 
#   filter(str_detect(gene_name, 'map'))
# 
# View(brain_f1_exp)

mapk_pathway = read_csv('mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  rename(gene_name = Symbol)

brain_div_plast_mapk = brain_f1_exp %>% 
  filter(gene_name %in% mapk_pathway$gene_name)

liver_div_plast_mapk = liver_f1_exp %>% 
  filter(gene_name %in% mapk_pathway$gene_name)

liver_div_plast_mapk %>% 
  select(gene_name) %>% 
  write_tsv('Liver_mapk_ecotemp_divergence.txt')




# skr hybrid 12 degree trans exp -------------------------------------------------------------

Gene_ID_12 = read_csv("~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_Transgressive_expression_12degrees.csv")

Trans_amb_hyb_12 = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_amb_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  # mutate(status = 'Outlier') %>% 
  inner_join(., 
             Gene_ID_12) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID) %>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

brain_amb_hyb_12_mapk = Trans_amb_hyb_12 %>% 
  filter(GeneID %in% mapk_pathway$gene_name)


Trans_geo_hyb_12 = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_geo_hyb_12_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  inner_join(., 
             Gene_ID_12) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID)%>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

brain_geo_hyb_12_mapk = Trans_geo_hyb_12 %>% 
  filter(GeneID %in% mapk_pathway$gene_name)


trans_mapk_overlap = inner_join(brain_amb_hyb_12_mapk, 
           brain_geo_hyb_12_mapk, 
           by = 'GeneID') %>% 
  select(GeneID)




# mapk trans exp 18 hybrids -----------------------------------------------


Gene_ID_18 = read_csv("~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_Transgressive_expression_18degrees.csv")

Trans_amb_hyb_18 = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_amb_hyb_18_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  # mutate(status = 'Outlier') %>% 
  inner_join(., 
             Gene_ID_18) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID) %>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

brain_amb_hyb_18_mapk = Trans_amb_hyb_18 %>% 
  filter(GeneID %in% mapk_pathway$gene_name)


Trans_geo_hyb_18 = read_csv('~/Parsons_Postdoc/SKR_Hybrid_Gene_expression/Brain_geo_hyb_18_div.csv')%>% 
  filter(adj.P.Val <= 0.05) %>% 
  inner_join(., 
             Gene_ID_18) %>% 
  dplyr::select(GeneID, 
                logFC, 
                adj.P.Val) %>% 
  rename(gene_ensembl = GeneID)%>% 
  inner_join(., 
             annotation_data) %>% 
  rename(GeneID = gene_name)

brain_geo_hyb_18_mapk = Trans_geo_hyb_18 %>% 
  filter(GeneID %in% mapk_pathway$gene_name)


trans_mapk_overlap_18 = inner_join(brain_amb_hyb_18_mapk, 
                                brain_geo_hyb_18_mapk, 
                                by = 'GeneID') %>% 
  select(GeneID)



# Hybrid divergence across temps MAPK -------------------------------------

mapk_disregulation_hybs = inner_join(trans_mapk_overlap, 
           trans_mapk_overlap_18, 
           by = 'GeneID')

