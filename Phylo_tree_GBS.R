######################################
## reimagining a phylogenetic tree
## 
## 02.20.2025
##
## Matt Brachmann (MKB) @phdmattyb
######################################

setwd('~/Parsons_Postdoc/Stickleback_Genomic/GBS_NJ_tree/')

library(tidyverse)
library(BiocManager)
library(ggtree)
library(ape)
library(tidytree)
library(ggtreeExtra)
library(ggnewscale)
library(ggstar)
library(reshape2)

tree_data = read.tree("GBS-NJ-TREEboot_main.nwk")

other_data = read.tree('GBS-NJ-TREEboot_other.nwk')
raxml_data = read.tree('GBS-2.raxml.support')

tree_data_tidy = as_tibble(tree_data) %>% 
  separate(col = label, 
           into = c('popeco', 
                    'ind'), 
           sep = '_') %>% 
  mutate(.data = .,
         Ecotype = as.factor(case_when(
           popeco == 'ASHNC' ~ 'Ambient',
           popeco == 'ASHNW' ~ 'Geothermal',
           popeco == 'BARN' ~ 'Geothermal',
           popeco == 'CSWY' ~ 'Ambient',
           popeco == 'GTS' ~ 'Geothermal',
           popeco == 'HERD' ~ 'Marine',
           popeco == 'HGRNS' ~ 'Ambient',
           popeco == 'KLFC' ~ 'Ambient',
           popeco == 'LITA' ~ 'Geothermal',
           popeco == 'LITP' ~ 'Geothermal',
           popeco == 'MYVC' ~ 'Ambient',
           popeco == 'MYVW' ~ 'Geothermal',
           popeco == 'NH' ~ 'Geothermal',
           popeco == 'NYPS' ~ 'Marine',
           popeco == 'OPNUR' ~ 'Geothermal',
           popeco == 'RKLTC' ~ 'Ambient',
           popeco == 'RKLTW' ~ 'Geothermal',
           popeco == 'RKR' ~ 'Geothermal',
           popeco == 'RKRC' ~ 'Ambient',
           popeco == 'SKAL' ~ 'Ambient',
           popeco == 'SKRC' ~ 'Ambient',
           popeco == 'SKRW' ~ 'Geothermal',
           popeco == 'STNST' ~ 'Geothermal',
           popeco == 'THNGC' ~ 'Ambient',
           popeco == 'THNGW' ~ 'Geothermal'))) %>% 
  mutate(.data = .,
         Population = as.factor(case_when(
           popeco == 'ASHNC' ~ 'ASHN',
           popeco == 'ASHNW' ~ 'ASHN',
           popeco == 'BARN' ~ 'BARN',
           popeco == 'CSWY' ~ 'GAR',
           popeco == 'GTS' ~ 'GTS',
           popeco == 'HERD' ~ 'HERD',
           popeco == 'HGRNS' ~ 'HGRNS',
           popeco == 'KLFC' ~ 'KLFC',
           popeco == 'LITA' ~ 'LITA',
           popeco == 'LITP' ~ 'LITP',
           popeco == 'MYVC' ~ 'MYV',
           popeco == 'MYVW' ~ 'MYV',
           popeco == 'NH' ~ 'NH',
           popeco == 'NYPS' ~ 'NYPS',
           popeco == 'OPNUR' ~ 'OPNUR',
           popeco == 'RKLTC' ~ 'RKLT',
           popeco == 'RKLTW' ~ 'RKLT',
           popeco == 'RKR' ~ 'RKR',
           popeco == 'RKRC' ~ 'RKR',
           popeco == 'SKAL' ~ 'SKAL',
           popeco == 'SKRC' ~ 'SKR',
           popeco == 'SKRW' ~ 'SKR',
           popeco == 'STNST' ~ 'STNST',
           popeco == 'THNGC' ~ 'THNG',
           popeco == 'THNGW' ~ 'THNG')))%>% 
  mutate(.data = .,
         Type = as.factor(case_when(
           popeco == 'ASHNC' ~ 'Sympatric',
           popeco == 'ASHNW' ~ 'Sympatric',
           popeco == 'BARN' ~ 'Sympatric',
           popeco == 'CSWY' ~ 'Allopatric',
           popeco == 'GTS' ~ 'Allopatric',
           popeco == 'HERD' ~ 'Marine',
           popeco == 'HGRNS' ~ 'Allopatric',
           popeco == 'KLFC' ~ 'Allopatric',
           popeco == 'LITA' ~ 'Sympatric',
           popeco == 'LITP' ~ 'Sympatric',
           popeco == 'MYVC' ~ 'Sympatric',
           popeco == 'MYVW' ~ 'Sympatric',
           popeco == 'NH' ~ 'Sympatric',
           popeco == 'NYPS' ~ 'Marine',
           popeco == 'OPNUR' ~ 'Allopatric',
           popeco == 'RKLTC' ~ 'Allopatric',
           popeco == 'RKLTW' ~ 'Allopatric',
           popeco == 'RKR' ~ 'Sympatric',
           popeco == 'RKRC' ~ 'Sympatric',
           popeco == 'SKAL' ~ 'Sympatric',
           popeco == 'SKRC' ~ 'Allopatric',
           popeco == 'SKRW' ~ 'Allopatric',
           popeco == 'STNST' ~ 'Sympatric',
           popeco == 'THNGC' ~ 'Sympatric',
           popeco == 'THNGW' ~ 'Sympatric'))) %>% 
  mutate(.data = .,
         LatLong = as.factor(case_when(
           popeco == 'ASHNC' ~ 'Sauderkroker',
           popeco == 'ASHNW' ~ 'Sauderkroker',
           popeco == 'BARN' ~ 'North-East',
           popeco == 'CSWY' ~ 'Sauderkroker',
           popeco == 'GTS' ~ 'North',
           popeco == 'HERD' ~ 'Marine',
           popeco == 'HGRNS' ~ 'North',
           popeco == 'KLFC' ~ 'South-West',
           popeco == 'LITA' ~ 'North-East',
           popeco == 'LITP' ~ 'North-East',
           popeco == 'MYVC' ~ 'North-East',
           popeco == 'MYVW' ~ 'North-East',
           popeco == 'NH' ~ 'North',
           popeco == 'NYPS' ~ 'Marine',
           popeco == 'OPNUR' ~ 'South-West',
           popeco == 'RKLTC' ~ 'West',
           popeco == 'RKLTW' ~ 'West',
           popeco == 'RKR' ~ 'North',
           popeco == 'RKRC' ~ 'North',
           popeco == 'SKAL' ~ 'North-East',
           popeco == 'SKRC' ~ 'Sauderkroker',
           popeco == 'SKRW' ~ 'Sauderkroker',
           popeco == 'STNST' ~ 'North',
           popeco == 'THNGC' ~ 'South-West',
           popeco == 'THNGW' ~ 'South-West'))) %>% 
  mutate(.data = .,
         Region = as.factor(case_when(
           popeco == 'ASHNC' ~ 'Lake 1',
           popeco == 'ASHNW' ~ 'Lake 1',
           popeco == 'BARN' ~ 'Lake 3',
           popeco == 'CSWY' ~ 'Lake 8',
           popeco == 'GTS' ~ 'Lake 9',
           popeco == 'HERD' ~ 'Marine',
           popeco == 'HGRNS' ~ 'Lake 10',
           popeco == 'KLFC' ~ 'Lake 12',
           popeco == 'LITA' ~ 'Lake 4',
           popeco == 'LITP' ~ 'Lake 4',
           popeco == 'MYVC' ~ 'Lake 2',
           popeco == 'MYVW' ~ 'Lake 2',
           popeco == 'NH' ~ 'Lake 3',
           popeco == 'NYPS' ~ 'Marine',
           popeco == 'OPNUR' ~ 'Lake 11',
           popeco == 'RKLTC' ~ 'Lake 7',
           popeco == 'RKLTW' ~ 'Lake 7',
           popeco == 'RKR' ~ 'Lake 3',
           popeco == 'RKRC' ~ 'Lake 3',
           popeco == 'SKAL' ~ 'Lake 4',
           popeco == 'SKRC' ~ 'Lake 6',
           popeco == 'SKRW' ~ 'Lake 6',
           popeco == 'STNST' ~ 'Lake 3',
           popeco == 'THNGC' ~ 'Lake 5',
           popeco == 'THNGW' ~ 'Lake 5'))) %>% 
  unite(sample_id, 
        c('popeco', 
          'ind'), 
        sep = '_', 
        remove = F) %>% 
  unite(popeco_ind, 
        c('popeco', 
          'ind'), 
        sep = '_', 
        remove = F)

# tree_data_tidy %>% 
#   group_by(popeco) %>% 
#   summarize(n = n()) %>% View()

phylo_data = as.treedata(tree_data_tidy)

# ggtree(phylo_data, 
#        layout="equal_angle")
# 
# ggtree(phylo_data, 
#        layout="daylight")

# ggtree(phylo_data, 
#        layout="daylight", 
#        branch.length = 'none')+
#   # geom_tiplab(aes(color = popeco), 
#   #             as_ylab=TRUE)
#   geom_point(aes(color = Population))
# +
#   geom_cladelab(node=489, 
#                 label="Marine", 
#                 angle=0, 
#                 fontsize=4) 
#   geom_text2(aes(label = popeco, 
#                  color = popeco), 
#              angle = 90,
#              size = 2)

# geom_tippoint(aes(colour = Population),
#                             alpha=0) +
# geom_tiplab(aes(colour = Population),
#             align=TRUE,
#             # linetype=3,
#             size=2,
#             linesize=0.2)+


population_pal2 = c('#277da1',
                   '#adb5bd',
                   '#9d4edd',
                   '#EDAE49',
                   '#b5e48c',
                   '#adb5bd',
                   '#adb5bd',
                   '#adb5bd',
                   '#adb5bd',
                   '#D1495B',
                   '#adb5bd',
                   '#b5e48c',
                   '#adb5bd',
                   '#adb5bd',
                   '#adb5bd',
                   '#adb5bd',
                   '#7ae7c7',
                   '#adb5bd',
                   '#adb5bd')

Geography_pal = c('#001219', 
          '#94d2bd', 
          '#0a9396',
         "#e9d8a6",
         "#9b2226",
         "#bb3e03")

region_pal = c('#277da1',
               '#3a5a40',
               '#344e41',
               '#606c38',
               '#D1495B',
               '#a2d2ff',
               '#dad7cd', 
               '#a3b18a',
               '#7ae7c7',
               '#588157',
               # '#ffafcc',
               '#9d4edd',
               '#EDAE49',
               '#b5e48c')

GBS_Tree = ggtree(phylo_data, 
       layout='circular', 
       aes(colour = Population)) + 
  xlim(-10, NA)+
  scale_color_manual(values = population_pal2)+
  geom_fruit(geom=geom_tile,
    mapping = aes(x = Ecotype, 
                y = ind, 
                fill = Ecotype), 
    color="black",) +
  scale_fill_manual(
    name="Ecotype",
    values=c("#003049", "#c1121f", "#b5e48c"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3)) +
  new_scale_fill()+
  geom_fruit(geom=geom_tile,
             mapping = aes(x = Region, 
                           y = ind, 
                           fill = Region), 
             color="black",) +
  scale_fill_manual(
    name="Geography",
   values = region_pal,
    na.translate=F,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3)) +
  # geom_fruit(geom=geom_point,
  #            mapping=aes(x=Type,
  #                        y=ind,
  #                        shape = Type),
  #            size=1,
  #            # starstroke=0,
  #            pwidth=0.1)+
  # new_scale_fill()+
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5.5),
    legend.spacing.y = unit(0.02, "cm"), 
    
  )+
  theme(legend.position = 'none')

GBS_Tree
ggsave('GBS_Phylo_tree_MKB_NoLegend_LakeRegion.tiff', 
       plot = GBS_Tree, 
       dpi = 'retina', 
       units = 'cm', 
       width = 50, 
       height = 50)



# PCA axes and adding to phylo tree ---------------------------------------

PCA_cov = read_table2('GBSfinal_unlinked-PCANGSD.cov',
         col_names = F) 


#We will also add a column with population assignments
pop <- c(rep("ASNHC",16),rep("ASNHW",16),rep("BARN",15),rep("CSWY",15),rep("GTS",14),rep("HERD",15),rep("HGRNS",13),rep("KLFC",12),rep("LITA",15),rep("LITP",15),rep("MYVC",16),rep("MYVW",16),rep("NH",14),rep("NYPS",15),rep("OPNUR",12),rep("RKLTC",15),rep("RKLTW",15),rep("RKR",13),rep("RKRC",11),rep("SKAL",15),rep("SKRC",14),rep("SKRW",16),rep("STNST",15),rep("THNGC",15),rep("THNGW",14))

#perform the pca using the eigen function.
e <- eigen(PCA_cov)  

#extract eigenvectors
eigenvectors = e$vectors 

#combine with the population assignments and plot
pca.vectors = as_tibble(cbind(pop, data.frame(eigenvectors))) 
pca = ggplot(data = pca.vectors, aes(x=X1, y=X2, colour = pop)) + geom_point()

#save
ggsave(filename = paste0(basedir, "Stats/PCA/pca2D_LDpruned_plot.pdf"), plot = pca) 



# admixture plot ----------------------------------------------------------

admix_pal = c('#3a5a40', ##lake 10
              '#277da1',#K1 ASHN
              '#9d4edd', ##CSWY
              '#D1495B', ## MYV
              '#344e41',
              '#b5e48c', ##Marine (HERD)
              'black',
              '#dad7cd', ## lake 4
              '#EDAE49', #GTS
              
              '#588157', ##lake 7
              '#a3b18a', ## SKR
              '#a2d2ff'# BARN group
) ##OTHER?
              

admix_meta_data = tree_data_tidy %>% 
  select(sample_id:Region) %>% 
  na.omit() %>% 
  arrange(popeco)


#The best K=11 run is the run 4
admix_data_init = read_table2('allindsK11run4.qopt', 
            col_names = F) %>% 
  rename(K1 = 1, 
         K2 = 2, 
         K3 = 3, 
         K4 = 4, 
         K5 = 5, 
         K6 = 6, 
         K7 = 7, 
         K8 = 8, 
         K9 = 9, 
         K10 = 10, 
         K11 = 11, 
         other = 12) %>% 
  dplyr::select(1:11)

populations<-c(rep("ASNHC",16),rep("ASNHW",16),rep("BARN",15),rep("CSWY",15),rep("GTS",14),rep("HERD",15),rep("HGRNS",13),rep("KLFC",12),rep("LITA",15),rep("LITP",15),rep("MYVC",16),rep("MYVW",16),rep("NH",14),rep("NYPS",15),rep("OPNUR",12),rep("RKLTC",15),rep("RKLTW",15),rep("RKR",13),rep("RKRC",11),rep("SKAL",15),rep("SKRC",14),rep("SKRW",16),rep("STNST",15),rep("THNGC",15),rep("THNGW",14))

admix_data = populations %>% 
  as_tibble() %>% 
  rename(popeco = value) %>% 
  bind_cols(.,
            admix_meta_data) %>% 
  dplyr::select(-popeco...1) %>% 
  rename(popeco = popeco...4) %>% 
  bind_cols(., 
            admix_data_init) %>% 
  dplyr::select(popeco, 
                ind, 
                K1:K11) %>% 
  unite(popeco_ind, 
        c('popeco', 
          'ind'), 
        sep = "_", 
        remove =F) %>% 
  dplyr::select(popeco_ind, K1:K11)

admix_data_reshaped = reshape2::melt(admix_data, 
     # id.vars = c('popeco', 
     #             'ind')) %>%
     id.vars = c('popeco_ind')) %>% 
  as_tibble() %>% 
  separate(popeco_ind, 
           into = c('popeco', 
                    'ind'), 
           remove = F) %>% 
  # unite(popeco_ind, 
  #       c('popeco', 
  #         'ind'), 
  #       sep = "_", 
  #       remove =F) %>% 
  rename(THEPOPS = popeco, 
         THEINDS = ind, 
         THECOMBO = popeco_ind)


admixture_plot = ggplot(data = admix_data_reshaped, 
       aes(x = reorder(THECOMBO, THEPOPS),
           y = value, 
           fill = variable, 
           group = THEPOPS))+
  geom_bar(stat = "identity", 
           width = 1)+
  scale_fill_manual(values = admix_pal)+
  # scale_fill_manual(values = magma(n = 4))+
  labs(x = 'Individuals', 
       y = 'Ancestry proportion')+
  theme(axis.text.y = element_text(color = 'black', 
                                   size = 12),
        axis.title.y = element_text(size = 14),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        ## can add xaxis labels if needed
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = -0.09,
                                   size = 6,
                                   color = 'black'),
        legend.position = 'none'
        )+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))

ggsave('GBS_admixture_plot.svg', 
       plot = admixture_plot, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 10)

# combine tree with admix plot --------------------------------------------

phylo_data = as.treedata(tree_data_tidy)

tibs_df_phylo = as_tibble(phylo_data) %>% 
  # rename(popeco_ind = sample_id)
  full_join(.,
            admix_data, 
            by = c('popeco_ind', 
                   'popeco', 
                   'ind')) %>% 
  dplyr::select(popeco_ind, 
                # popeco, 
                # ind, 
                K1:K11)

admix_plot_data = reshape2::melt(tibs_df_phylo, 
               id.vars = c('popeco_ind')) %>% 
  as_tibble() %>% 
  separate(col = popeco_ind, 
           c('popeco', 
             'ind'), 
           sep = '_', 
           remove = F)
# %>% 
  # unite(popeco_ind, 
  #       c('popeco', 
  #         'ind'), 
  #       sep = "_", 
  #       remove =F) %>% 
  # rename(THEPOPS = popeco, 
  #        THEINDS = ind, 
  #        THECOMBO = popeco_ind)
  # 



GBS_Tree+
  # coord_flip()+
  new_scale_fill()+
  geom_fruit(
    data = admix_plot_data,
    geom=geom_bar,
             mapping = aes(x = value,
               # x = reorder(THECOMBO, THEPOPS),
                           y = ind, 
                           fill = variable), 
    stat = "identity",
    orientation = "y") +
  scale_fill_manual(
    name="Admixture",
    values=col_vector,
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3))



# Phylo tree WGS data -----------------------------------------------------

library(vcfR)
library(phytools)
library(fastreeR)
library(dartRverse)

vcf = read.vcfR('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/stickleback_filtered_vcf.vcf')
vcf_names = colnames(vcf@gt) %>% 
  as_tibble() %>% 
  slice(-1)


metadata = read_csv('~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/stickleback_identifiers.csv') %>% 
  bind_cols(., 
            vcf_names) %>% 
  rename(vcf_names = value)
# %>%
#   mutate(individual_id2 = individual_id) %>% 
#   unite(col = individual_id, 
#         c('individual_id', 
#           'individual_id2'), 
#         sep = '_')
genlight_data <- vcfR2genlight(vcf)

# dist_matrix = vcf2dist(inputFile = '~/Parsons_Postdoc/Stickleback_Genomic/vcf_filter/stickleback_filtered_vcf.vcf',
#          outputFile = 'stickleback_dist_matrix.txt')

# library(poppr)
# 
# nei.dist(genlight_data, warning = TRUE)

WGS_dist_mat = dist(genlight_data)

WGS_tree = upgma(WGS_dist_mat)

WGS_tree_data = as_tibble(WGS_tree) %>%
  dplyr::rename(vcf_names = label) %>% 
  full_join(., 
             metadata, 
             by = 'vcf_names') %>% 
  mutate(.data = .,
         Ecotype = as.factor(case_when(
           population == 'ASHNC' ~ 'Ambient',
           population == 'ASHNW' ~ 'Geothermal',
           population == 'CSWY' ~ 'Ambient',
           population == 'GTS' ~ 'Geothermal',
           population == 'MYVC' ~ 'Ambient',
           population == 'MYVW' ~ 'Geothermal',
           population == 'SKRC' ~ 'Ambient',
           population == 'SKRW' ~ 'Geothermal'))) %>% 
  mutate(.data = .,
         Population2 = as.factor(case_when(
           population == 'ASHNC' ~ 'ASHN',
           population == 'ASHNW' ~ 'ASHN',
           population == 'CSWY' ~ 'GAR',
           population == 'GTS' ~ 'GTS',
           population == 'MYVC' ~ 'MYV',
           population == 'MYVW' ~ 'MYV',
           population == 'SKRC' ~ 'SKR',
           population == 'SKRW' ~ 'SKR')))

WGS_tree = as.treedata(WGS_tree_data)


WGS_pop_pal1 = c('#277da1',
                    '#9d4edd',
                    '#EDAE49',
                    '#D1495B',
                    '#7ae7c7', 
                    'black')

WGS_pop_tree = ggtree(WGS_tree, 
       layout="daylight", 
       branch.length = 'none', 
       size = 1)+
  geom_tippoint(size=3, 
                aes(colour = Population2))+
  scale_color_manual(values = WGS_pop_pal1)+
  theme(legend.position = 'none')

ggsave('WGS_Phylo_tree_MKB_NoLegend.svg', 
       plot = WGS_pop_tree, 
       dpi = 'retina', 
       units = 'cm', 
       width = 10, 
       height = 20)



WGS_pop_pal2 = c('#a2d2ff', 
                 '#ffafcc', 
                 '#9d4edd', 
                 '#EDAE49', 
                 '#003049', 
                 '#780000', 
                 '#2a9d8f', 
                 '#e76f51', 
                 'black')



ggtree(WGS_tree, 
       layout="daylight", 
              branch.length = 'none', 
       aes(colour = population))+
  geom_tippoint(size=3, 
                aes(colour = population))+
  scale_color_manual(values = WGS_pop_pal2)

  #   new_scale_fill()+
  # geom_fruit(geom = geom_tile,
  #            mapping = aes(x = Ecotype, 
  #                          y = vcf_names, 
  #                          fill = Ecotype), 
  #            color="black") +
  # scale_fill_manual(
  #   name="Ecotype",
  #   values=c("#003049", "#c1121f"),
  #   na.translate=FALSE,
  #   guide=guide_legend(keywidth=0.5,
  #                      keyheight=0.5,
  #                      order=3)) 

ggtree(WGS_tree, 
       layout='circular', 
       aes(colour = population)) + 
  xlim(-10, NA)
