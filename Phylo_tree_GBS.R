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

ggtree(phylo_data, 
       layout='circular', 
       aes(colour = Population)) + 
  xlim(-10, NA)+
  ggtreeExtra::geom_fruit(geom=geom_tile,
    mapping = aes(x = Ecotype, 
                y = ind, 
                fill = Ecotype)) +
  scale_fill_manual(
    name="Ecotype",
    values=c("#595959", "#B30000", "#020099", "#E6E6E6"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3)) +
  new_scale_fill()+
  geom_fruit(geom=geom_point,
    mapping=aes(x=Type, 
                y=ind, 
                fill=Type, 
                color = Type, 
                shape = Type),
    size=1,
    # starstroke=0,
    pwidth=0.1,
    grid.params=list(
      linetype=3,
      size=0.2))


  geom_point(aes(color = sample_id))
  geom_text2(aes(label = popeco, 
                 color = popeco), 
             angle = 90, 
             size = 2)
  # geom_text2(aes(label = popeco,
  #                color = popeco,
  #                subset = !is.na(as.numeric(ind))),
  #            vjust = 1.2,
  #            hjust = 1.2)+
  scaleClade(p, 23, .2) %>% 
  collapse(23, 'min', fill="darkgreen")  

  # geom_hilight(mapping=aes(subset = node %in% c(10, 12), 
  #                          fill = S),
  #              type = "gradient", 
  #              gradient.direction = 'rt',
  #              alpha = .8)
  # geom_label(aes(x=branch, 
  #                label=Ecotype, 
  #                fill = Ecotype))

  # geom_treescale(offset = -1)
  # theme_tree2()

ggtree(phylo_data, 
       ladderize = F)+
  geom_text2(aes(label = popeco,
                 subset = !is.na(as.numeric(ind))),
             vjust = 1.2,
             hjust = 1.2)


