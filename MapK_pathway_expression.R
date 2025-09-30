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

liver_mapk_genes = liver_div_plast_mapk %>% 
  select(ensemble_name,, 
         gene_name) %>% 
  rename(GeneID = ensemble_name)

# limma -------------------------------------------------------------------

liver_exp_data = read_table2('~/Parsons_Postdoc/Bethany_gene_expression/F1_lab_liver_GE_over10M.tsv') %>% 
  rename(ensemble_id = GeneID)

liver_meta_data = read_csv("~/Parsons_Postdoc/Bethany_gene_expression/F1_labfish_sampledata10M_liver.csv") %>% 
  unite(col = ecotemp, 
        c('popeco', 
          'temp'), 
        sep = '_', 
        remove = F) %>% 
  unite(col = ecotemp_simple, 
        c('ecotype', 
          'temp'), 
        sep = '_', 
        remove = F)

gene_name_data = annotation_data %>% 
  ungroup() %>% 
  select(ensemble_id, 
         gene_name)

liver_exp_data = inner_join(gene_name_data,
           liver_exp_data,
           by = 'ensemble_id') %>% 
  select(-ensemble_id)

liver_mapk_counts = liver_exp_data %>% 
  filter(gene_name %in% mapk_pathway$gene_name)


liver_mapk_dge = DGEList(liver_mapk_counts)
liver_mapk_norm = calcNormFactors(liver_mapk_dge)

mm = model.matrix(~0 + ecotemp_simple, 
                  data = liver_meta_data)

liver_mapk_keep = filterByExpr(liver_mapk_norm, 
                          min.count = 2,
                          mm)
sum(liver_mapk_keep) # number of genes retai

liver_mapk_keep = liver_mapk_norm[liver_mapk_keep,]

## EdgeR model
liver_mapk_dispersion = estimateDisp(liver_mapk_keep, 
                                mm) 

contrast = makeContrasts(cold_plast = ecotemp_simpleC_12oC - ecotemp_simpleC_18oC, 
                         warm_plast = ecotemp_simpleW_12oC - ecotemp_simpleW_18oC,
                         eco_div_12 = ecotemp_simpleC_12oC - ecotemp_simpleW_12oC, 
                         eco_div_18 = ecotemp_simpleC_18oC - ecotemp_simpleW_18oC, 
                         levels = mm)

# contrast = makeContrasts(ASHN_cold_plast = ecotempASHNC_12oC - ecotempASHNC_18oC, 
#                          ASHN_warm_plast = ecotempASHNW_12oC - ecotempASHNW_18oC,
#                          ASHN_eco_div_12 = ecotempASHNC_12oC - ecotempASHNW_12oC, 
#                          ASHN_eco_div_18 = ecotempASHNC_18oC - ecotempASHNW_18oC,
#                          MYV_cold_plast = ecotempMYVC_12oC - ecotempMYVC_18oC, 
#                          MYV_warm_plast = ecotempMYVW_12oC - ecotempMYVW_18oC, 
#                          MYV_eco_div_12 = ecotempMYVC_12oC - ecotempMYVW_12oC, 
#                          MYV_eco_div_18 = ecotempMYVC_18oC - ecotempMYVW_12oC,
#                          SKR_cold_plast = ecotempSKRC_12oC - ecotempSKRC_18oC, 
#                          SKR_warm_plast = ecotempSKRW_12oC - ecotempSKRW_18oC, 
#                          SKR_eco_div_12 = ecotempSKRC_12oC - ecotempSKRW_12oC, 
#                          SKR_eco_div_18 = ecotempSKRC_18oC - ecotempSKRW_18oC,
#                          levels = mm)

liver_mapk_glm_fit = glmQLFit(liver_mapk_dispersion, 
                         # contrast = ecotype.div.brain,
                         design = mm)

liver_mapk_glm_test = glmQLFTest(liver_mapk_glm_fit, 
                            contrast = contrast)

liver_mapk_edger_results = topTags(liver_mapk_glm_test, 
                              n = Inf,
                              adjust.method = 'bonferroni', 
                              p.value = 0.05)

liver_mapk_edger_results$table %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  write_csv('Liver_mapk_edger_results.csv')


# ## limma model
liver_mapk_voom = voom(liver_mapk_keep, 
                  mm, 
                  plot = T)
# 
liver_mapk_fit_limma <- limma::lmFit(liver_mapk_voom, 
                                contrast = contrast,
                                design=mm)
liver_mapk_fit_limma_contrast = contrasts.fit(liver_mapk_fit_limma, 
                                         contrasts = contrast)


liver_mapk_fit_ebayes = eBayes(liver_mapk_fit_limma_contrast)

liver_mapk_limma_results = topTable(liver_mapk_fit_ebayes, 
                               n = Inf, 
                               adjust.method = 'bonferroni', 
                               p.value = 0.05)


# DESEQ MAPK gene expression ----------------------------------------------

library(DESeq2)

# liver_exp_data = read_table2('~/Parsons_Postdoc/Bethany_gene_expression/F1_lab_liver_GE_over10M.tsv') %>% 
#   rename(ensemble_id = GeneID)

liver_mapk_counts = read_table2('Liver_mapk_count_data.txt') %>% 
  select(-gene_name)


liver_meta_data = read_csv("~/Parsons_Postdoc/Bethany_gene_expression/F1_labfish_sampledata10M_liver.csv") %>% 
  unite(col = ecotemp, 
        c('popeco', 
          'temp'), 
        sep = '_', 
        remove = F) %>% 
  unite(col = ecotemp_simple, 
        c('ecotype', 
          'temp'), 
        sep = '_', 
        remove = F)

# dds = DESeqDataSetFromMatrix(countData = liver_mapk_counts, 
#                              colData = liver_meta_data, 
#                              design = ~ ecotemp_simple)
# dds = estimateSizeFactors(dds)
# dds = estimateDispersions(dds)
# 
# vst = getVarianceStabilizedData(dds)
# 
# vst %>% 
#   as_tibble() %>% 
#   write_tsv('Liver_mapk_Normalized_counts.txt')

liver_mapk_counts_norm = read_tsv('Liver_mapk_Normalized_counts.txt')
liver_mapk_counts = read_table2('Liver_mapk_count_data.txt') %>% 
  select(-gene_name)
liver_mapk_genes = read_table2('Liver_mapk_count_data.txt') %>% 
  select(gene_name)



dds = model.matrix(~ ecotype + temp + ecotype*temp, 
                   liver_meta_data)
# colnames(dds)
# unname(dds)
# all.zero <- apply(dds, 2, function(x) all (x==0))
# all.zero
# idx <- which(all.zero)
# dds <- dds[,-idx]
# unname(dds)
# lapply(liver_meta_data, class)

##now model it, following suggestions for large models with rows which do not converge in beta##
object = DESeqDataSetFromMatrix(countData = liver_mapk_counts, 
                                colData = liver_meta_data, 
                                design = ~ ecotype + temp + ecotype*temp, 
                                ignoreRank = T)
##filter out rows where gene is expressed in less than half of the samples##
object <- estimateSizeFactors(object)
# nc <- counts(object, normalized=T)
# object <- object[ rowSums(nc > 0) >= 163 ]

##actually run the analysis; do it in two steps to increase power and get more rows to converge##
object <- estimateDispersions(object)
#save(colData, countData, dds, object, file = "Stickeblack-Frozed.RData")
#load("Stickeblack-Frozed.RData")
dds1 <- nbinomWaldTest(object, maxit=999)


res <- results(dds1)

res %>% 
  as_tibble() %>%
  bind_cols(., 
            liver_mapk_genes) %>% 
  filter(pvalue < 0.05)


##Look at Infection##
resinfect <- results(dds1,contrast=list(c("worm_presentTRUE")))
resorderd <-resinfect[order(resinfect$padj),] 
head(resorderd, 2370)
write.csv(resorderd, file = "DESeq_Infection_Full_StringParam_new.csv") ##2369 affected


# norm count graph --------------------------------------------------------

liver_mapk_counts_norm = read_tsv('Liver_mapk_Normalized_counts.txt')
liver_mapk_genes = read_table2('Liver_mapk_count_data.txt') %>% 
  select(gene_name)

liver_mapk_norm_counts_log = bind_cols(liver_mapk_genes, 
                                   liver_mapk_counts_norm) %>% 
  # filter(gene_name %in% c('braf', 
  #                         'prkcg')) %>% 
  mutate(across(where(is.numeric), ~ log(.x)))

liver_mapk_log_norm_counts_clean = liver_mapk_norm_counts_log %>% 
  pivot_longer(cols = !(gene_name), 
               names_to = 'samples', 
               values_to = 'log_count') %>% 
  filter(gene_name %in% c('braf', 
                          'prkcg')) %>% 
  separate(col = samples, 
           into = c('sample_num', 
                    'ecopop', 
                    'temp'), 
           sep = '_') %>% 
  select(-sample_num) %>% 
  separate(col = temp, 
           into = c('temp', 
                    'fam', 
                    'num'))%>% 
  mutate(.data = .,
         ecotype = as.factor(case_when(
           ecopop == 'ASHNW' ~ 'Warm',
           ecopop == 'ASHNC' ~ 'Cold',
           ecopop == 'SKRW' ~ 'Warm', 
           ecopop == 'SKRC' ~ 'Cold', 
           ecopop == 'MYVW' ~ 'Warm', 
           ecopop == 'MYVC' ~ 'Cold'))) %>% 
  mutate(.data = ., 
         pop = as.factor(case_when(
           ecopop == 'ASHNW' ~ 'ASHN', 
           ecopop == 'ASHNC' ~ 'ASHN', 
           ecopop == 'SKRW' ~ 'SKR', 
           ecopop == 'SKRC' ~ 'SKR', 
           ecopop == 'MYVW' ~ 'MYV', 
           ecopop == 'MYVC' ~ 'MYV'
         ))) %>% 
  unite(col = ecotemp, 
        c('ecotype', 
          'temp'), 
        sep = '_', 
        remove = F)

liver_mapk_means = liver_mapk_log_norm_counts_clean %>% 
  group_by(gene_name,
           temp, 
           ecotype) %>% 
  summarize(mean_expression = mean(log_count))
  
liver_mapk_log_norm_counts_clean %>% 
  ggplot(aes(x = temp, 
             y = log_count, 
             col = ecotype))+
  geom_point()+
  geom_point(data = liver_mapk_means,
             aes(x = temp, 
                 y = mean_expression), 
             col = 'black') + 
  facet_grid(~gene_name)


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
