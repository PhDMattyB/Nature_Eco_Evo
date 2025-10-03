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

col_vector<-c("black","#d3436e","snow","seagreen","#f8765c","darkgoldenrod4",
              "gray","#fcfdbf","yellowgreen","#3683d3","#5f187f","navy")

#The best K=11 run is the run 4
admix<-t(as.matrix(read.table("allindsK11run4.qopt")))
head(admix[, 1:10])

admix_data = read_table2('allindsK11run4.qopt', 
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
  rename(Population = value) %>% 
  bind_cols(., 
            admix_data)


admix_data_reshaped = reshape2::melt(admix_data, 
     id.vars = c('Population')) %>% 
  as_tibble()


ggplot(data = admix_data_reshaped, 
       aes(x = Population,
           y = value, 
           fill = variable, 
           group = Population))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = col_vector)+
  scale_x_discrete(guide = guide_axis(n.dodge = 5))+
  scale_y_continuous(expand = c(0,0))


barplot(admix, col=col_vector, space=0, border=NA, ylab="Admixture", xlab="Populations", main="Icelandic Sticklebacks (K=11)")

abline(v=c(16,32,47,62,76,91,104,116,131,146,162,178,192,207,219,234,249,262,273,288,302,318,333,348),lty=5,lwd=2, col="white")
#text(c(1,17,33,48,63,77,92,105,117,132,147,163,179,193,208,220,235,250,263,274,289,303,319,334,349),-0.05,unique(populations),xpd=T)

#relocate the position of the groups in the plot to follow the same order as placed in the Map
geolocate<-admix[,c(1:16,17:32,33:47,179:192,319:333,250:262,263:273,147:162,163:178,117:131,132:146,274:288,
                    334:348,349:362,289:302,303:318,48:62,92:104,63:76,220:234,235:249,208:219,105:116,77:91,193:207)]

barplot(geolocate, col=col_vector, space=0, border=NA, ylab="Admixture", xlab="Populations", main="Icelandic Sticklebacks (K=11)")

abline(v=c(16,32,47,61,76,89,100,116,132,147,162,177,192,206,220,236,251,264,278,293,308,320,332,347),lty=5,lwd=2, col="white")

#text(c(1,17,33,48,62,78,91,105,120,134,147,158,173,189,205,220,235,250,265,279,294,309,321,333,348),-0.05,unique(populations),xpd=T)
