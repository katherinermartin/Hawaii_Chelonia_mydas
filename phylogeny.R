# Import tree file for Bayesian phylogeny of Cape Verde (Stiebens et al 2013), Florida (Martin et al 2022) and Hawaii Adkins (2022) MHC class I alleles from C. carett and C. mydas

library(ggtree) # for plotting
library(treeio) # for read.beast function
library(ape)
library(ggplot2) # for plotting
library(phytools)
library(svglite) # to save
library(tidyverse)


tree <-read.beast("infile.nex.con.tre")

df <- read.csv("Cc_Cm_CV_FL_HI_allele_counts_by_species_locality_zeroes_ones.csv")

# make the haplotype column the name of the rows
df <- df %>% remove_rownames %>% column_to_rownames(var="allele")
head(df)


# turn the matrix of zeroes and ones into the respective species and location:

df <- df %>% mutate_at("Cc_FL", str_replace, "1", "Cc_FL")
df <- df %>% mutate_at("Cc_CV", str_replace, "1", "Cc_CV")
df <- df %>% mutate_at("Cm_FL", str_replace, "1", "Cm_FL")
df <- df %>% mutate_at("Cm_HI", str_replace, "1", "Cm_HI")

# make all zeroes NAs in the dataframe
df[df == 0] <- NA

species_locations <- df %>% select(c("Cc_FL", "Cc_CV", "Cm_FL", "Cm_HI"))


tree <- read.mrbayes("infile.nex.con.tre")
tree # this has 137 tips.

# drop the outgroup tips:
tree <- treeio::drop.tip(tree, c("KF032390_Gallus_gallus", "KF466478_Tympanuchus_cupido"))

# Initialize tree: tree aesthetics
p1 <- ggtree(tree, # tree read in
             ladderize = TRUE,
             right = TRUE, 
             size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)
p1

# add internal nodes colored by posterior probability
p2 <- p1 +
  geom_nodepoint(color = "black", fill = "black", aes(subset = prob >= 0.995), size = 1, shape = 21) + #posterior prob between 0.995- 1
  geom_nodepoint(color = "black", fill = "#808080", aes(subset = prob < 0.995 & prob >= 0.945), size = 1, shape = 21) + #posterior prob between 0.945-0.994
  geom_nodepoint(color = "black", fill = "#dcdcdc", aes(subset = prob < 0.945 & prob >= 0.895), size = 1, shape = 21) # posterior prob between 0.895-0.944
p2

p3 <- p2 + geom_tiplab(color="black", size = 1.5, hjust=-.2) # adds tip labels
p3 # tree with correctly colored nodes and labels too.


p4 <- gheatmap(p2, # previous tree to use, without labels
               species_locations, # variable to create a heatmap from
               width=0.20, # width of bars
               color = "NA", # color of bar borders
               colnames=FALSE, # no column names
               offset = 0) + #distance between heatmap and tips
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3)) +
  scale_fill_manual(
    values = c("#ff8c00",
               "#ffcc00",
               "#008000",
               "#00cdcd"),
    breaks = c("Cc_FL",
               "Cc_CV",
               "Cm_FL",
               "Cm_HI"),
    na.value = "white", name = "species") # have to put name of legend here
p4

ggsave("", p4) # remove root, add legend, add supertype labels


# Maximum likelihood phylogeny
### load in  IQTree maximum likelihood phylogeny from https://www.hiv.lanl.gov/cgi-bin/IQTREE/iqtree.cgi; potential supplemental figure

tree <- read.tree("Cm_Hawaii_IQTree.nwk")

p <- ggtree(tree, ladderize= TRUE, right = TRUE) +
  geom_tiplab(color="black", size = 2, hjust = -0.2) +
  geom_nodelab(aes(label = label), size = 3, hjust = 1)
p # all bootstrap values labeled


p <- ggtree(tree, ladderize= TRUE, right = TRUE) +
  geom_tiplab(color="black", size = 1.8, hjust = -0.2) +
  geom_label2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) >= 90))
p # bootstrap values above 90 but in large, wonky labels


p <- ggtree(tree, ladderize= TRUE, right = TRUE) +
  geom_tiplab(color="black", size = 2, hjust = -0.2) +
  geom_nodepoint(color = "black", fill = "black", aes(subset = label >= 90), size = 1, shape = 21)
p # bootstrap labels above 90 labeled with solid black dot


ggsave("", p) # remove root, add legend, add supertype labels
