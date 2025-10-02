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

tree_data = read.tree("GBS-NJ-TREEboot_main.nwk")

other_data = read.tree('GBS-NJ-TREEboot_other.nwk')


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
           popeco == 'THNGW' ~ 'Sympatric')))

# tree_data_tidy %>% 
#   group_by(popeco) %>% 
#   summarize(n = n()) %>% View()

phylo_data = as.treedata(tree_data_tidy)


ggtree(phylo_data, 
       ladderize = T)+
  geom_text2(aes(label = popeco,
                 subset = !is.na(as.numeric(ind))),
             vjust = 1.2,
             hjust = 1.2)


