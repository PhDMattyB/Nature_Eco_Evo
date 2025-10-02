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
  unite(sample_id, 
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

GBS_Tree = ggtree(phylo_data, 
       layout='circular', 
       aes(colour = Population)) + 
  xlim(-10, NA)+
  scale_color_manual(values = population_pal2)+
  geom_fruit(geom=geom_tile,
    mapping = aes(x = Ecotype, 
                y = ind, 
                fill = Ecotype)) +
  scale_fill_manual(
    name="Ecotype",
    values=c("#003049", "#c1121f", "#b5e48c"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3)) +
  new_scale_fill()+
  geom_fruit(geom=geom_point,
    mapping=aes(x=Type, 
                y=ind, 
                shape = Type),
    size=1,
    # starstroke=0,
    pwidth=0.1)+
  theme(legend.position = 'none')


ggsave('GBS_Phylo_tree_MKB_sanslegend.svg', 
       plot = GBS_Tree, 
       dpi = 'retina', 
       units = 'cm', 
       width = 30, 
       height = 30)
