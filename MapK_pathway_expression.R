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


# annotation data ---------------------------------------------------------


gene_annotation = read_tsv('stickleback_v5_ensembl_genes.gff3.gz', 
                           col_names = F, 
                           skip = 1) %>% 
  # filter(X3 %in% c('gene', 
  #                  'exon', 
  #                  'CDS')) %>% 
  group_by(X1) %>% 
  arrange(X4, 
          X5) %>% 
  ## arrange each gene by its start and end points on each chromosome
  mutate(mid = X4 + (X5-X4)/2) %>% 
  dplyr::select(X1, 
                X3:X5, 
                X9:mid) %>% 
  rename(chromosome = X1, 
         feature = X3, 
         start = X4, 
         end = X5, 
         gene_id = X9, 
         position = mid) %>% 
  na.omit()

gene_annotation %>% 
  arrange(chromosome) %>%
  # filter(feature == 'gene')
  group_by(feature) %>% 
  distinct(feature)

annotation_data = gene_annotation %>% 
  # filter(chromosome == 'chrXXI') %>% 
  separate(col = gene_id, 
           into = c('ensemble_id', 
                    'gene_name', 
                    'parent_code'), 
           sep = ';') %>% 
  separate(col = ensemble_id, 
           into = c('garbage', 
                    'ensemble_id'), 
           sep = '=') %>% 
  dplyr::select(-garbage) %>% 
  separate(col = gene_name, 
           into = c('Garbage', 
                    'gene_name'), 
           sep = '=') %>% 
  dplyr::select(-Garbage) %>% 
  separate(col = parent_code, 
           into = c('garbage', 
                    'parent_gene_name'), 
           sep = '=') %>% 
  dplyr::select(-garbage) %>% 
  filter(feature == 'gene') 

# %>% 
#   ungroup(chromosome) %>% 
#   dplyr::select(gene_name) %>% 
#   distinct(gene_name) %>% 
#   filter(!if_any(everything(), ~ grepl('ENSG', .))) %>% 
#   mutate_all(as.character) %>% 
#   as.data.frame()



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
  rename(ensemble_id = GeneID) %>% 
  # rename(gene_ensembl = GeneID) %>% 
  inner_join(., 
             annotation_data, 
             by = 'ensemble_id') %>% 
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
  rename(ensemble_id = GeneID) %>% 
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
  rename(ensemble_id = GeneID) %>% 
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
  rename(ensemble_id = GeneID) %>% 
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



# Overlap 2 rna-seq datasets ----------------------------------------------
liver_div_plast_mapk %>% 
  dplyr::select(gene_name) %>% 
  rename(GeneID = gene_name) %>% 
  inner_join(mapk_disregulation_hybs)

liver_div_plast_mapk %>% 
  dplyr::select(gene_name) %>% 
  rename(GeneID = gene_name) %>% 
  inner_join(trans_mapk_overlap)

liver_div_plast_mapk %>% 
  dplyr::select(gene_name) %>% 
  rename(GeneID = gene_name) %>% 
  inner_join(trans_mapk_overlap_18)

# Liver plast hybrid expression plot--------------------------------------------------

inner_join(brain_amb_hyb_18_mapk, 
           brain_geo_hyb_18_mapk, 
           by = 'GeneID') %>% 
  filter(GeneID %in% c('braf', 
                           'prkcg')) 

amb_hyb_18_select_genes = brain_amb_hyb_18_mapk %>% 
  filter(GeneID %in% c('braf', 
                       'prkcg')) %>% 
  mutate(temp = '18', 
         div = 'Ambient-Hybrid')

geo_hyb_18_select_genes = brain_geo_hyb_18_mapk %>% 
  filter(GeneID %in% c('braf', 
                       'prkcg')) %>% 
  mutate(temp = '18', 
         div = 'Geothermal-Hybrid')

hybrid_div_18 = bind_rows(amb_hyb_18_select_genes, 
          geo_hyb_18_select_genes) %>% 
  unite(col = 'div_temp', 
        c('div', 
          'temp'), 
        sep = '_', 
        remove = F)


amb_hyb_12_select_genes = brain_amb_hyb_12_mapk %>% 
  filter(GeneID %in% c('braf', 
                       'prkcg')) %>% 
  mutate(temp = '12', 
         div = 'Ambient-Hybrid')

geo_hyb_12_select_genes = brain_geo_hyb_12_mapk %>% 
  filter(GeneID %in% c('braf', 
                       'prkcg')) %>% 
  mutate(temp = '12', 
         div = 'Geothermal-Hybrid')

hybrid_div_12 = bind_rows(amb_hyb_12_select_genes, 
                          geo_hyb_12_select_genes) %>% 
  unite(col = 'div_temp', 
        c('div', 
          'temp'), 
        sep = '_', 
        remove = F)

hybrid_div_key_genes = bind_rows(hybrid_div_12, 
                                 hybrid_div_18)

hyb_expression_cols = c('#003049', 
                        '#c1121f')

mapk_hyb_gene_expression = hybrid_div_key_genes %>% 
  mutate_if(is.character, toupper) %>%
  ggplot(aes(x = temp, 
             y = logFC, 
             col = div, 
             group = div))+
  geom_point(size = 3)+
  geom_line()+
  facet_grid(~GeneID)+
  scale_color_manual(values = hyb_expression_cols)+
  theme(legend.position = 'none', 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 14), 
        axis.text = element_text(size = 12), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        # panel.grid.major = element_blank(), 
        strip.background = element_rect(fill = 'white'))

ggsave('Mapk_expression_hybrid_genes.svg', 
       plot = mapk_hyb_gene_expression, 
       dpi = 'retina', 
       units = 'cm', 
       width = 10, 
       height = 10)

# Liver GLMER GO results --------------------------------------------------

liver_glmer_GO = read_csv('Liver_Mapk_GO.csv') %>% 
  dplyr::select(GO, 
                `Z-score`, 
                LogP, 
                Enrichment, 
                `Log(q-value)`, 
                Category, 
                Description) %>% 
  unite(GO_cat, 
        c('GO', 
          'Description'), 
        sep = ' - ', 
        remove = F)

db_col_pal = c(
  # '#031d44', 
               '#04395e', 
               '#70a288',
               '#dab785', 
               '#d5896f')

liver_glmer_go_plot = liver_glmer_GO %>% 
  # mutate_if(is.character, toupper) %>% 
  arrange(`Z-score`) %>% 
  # mutate(Category = factor(Category,levels = Category)) %>%
  mutate(GO_cat = factor(GO_cat, unique(GO_cat))) %>% 
  # rename(Database = database) %>% 
  ggplot(aes(y = GO_cat,
             x = `Z-score`, 
             fill = Category), 
         col = 'black')+
  # ggplot(aes(y = Term, 
  #            x = log_adj_pval, 
  #            col = Combined.Score, 
  #            fill = Combined.Score))+
  geom_col()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  scale_fill_manual(values = db_col_pal)+
  labs(x = 'Log q-value')+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size = 12), 
        legend.position = 'none')

ggsave('Metascape_GO_Term_Significant_plot_mapk_signalling.svg', 
       plot = liver_glmer_go_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 15)

# FST outliers MAPK pathway -----------------------------------------------

# ASHN_outliers = read_csv('ASHN_FST_0.5%_outlier_genes.csv')
# MYV_outliers = read_csv('MYV_FST_0.5%_outlier_genes.csv')
# SKR_outliers = read_csv('SKR_FST_0.5%_outlier_genes.csv')
# GTS_CSWY_outliers = read_csv('GTS_CSWY_FST_0.5%_outlier_genes.csv')

WC_outliers = read_csv("WC_FST_0.5%_outlier_genes.csv")

ASHN_outliers2 = read_csv('ASHN_FST_0.5%_outlier_genes_NEW.csv')
MYV_outliers2 = read_csv('MYV_FST_0.5%_outlier_genes_NEW.csv')
SKR_outliers2 = read_csv('SKR_FST_0.5%_outlier_genes_NEW.csv')
GTS_CSWY_outliers2 = read_csv('GTS_CSWY_FST_0.5%_outlier_genes_NEW.csv')



mapk_pathway = read_csv('mapk_pathway_genelist_fixed.csv')%>%
  mutate_all(., .funs = tolower) %>% 
  rename(gene_name = Symbol)


ASHN_mapk = inner_join(ASHN_outliers2, 
           mapk_pathway, 
           by = 'gene_name')

MYV_mapk = inner_join(MYV_outliers2, 
                      mapk_pathway, 
                      by = 'gene_name')

SKR_mapk = inner_join(SKR_outliers2, 
                      mapk_pathway, 
                      by = 'gene_name')

GTS_CSWY_mapk = inner_join(GTS_CSWY_outliers2, 
                           mapk_pathway, 
                           by = 'gene_name')

WC_mapk = inner_join(WC_outliers, 
                     mapk_pathway, 
                     by = 'gene_name')
