library(ggplot2)
library(dplyr)
#devtools::install_github("tidyverse/multidplyr")
library(multidplyr)
library(forcats)


lca = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  lastuni = tail(names(countshnames[countshnames==1]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni)
  return(ret)
}

# read file 
blasthits = read.table('blasthits2.tab')


# set up a new cluster of R sessions
cluster <- new_cluster(8)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca')


#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca(.)}) %>% # Use do to apply function across groups
  collect() %>% #bring it all back together
  slice(1) %>% #keep only first row in each group (one per read)
  group_by(last_common) %>% # group by LCA taxon
  summarize(lca_count = n()) %>% # count num reads (rows) for each LCA
  arrange(desc(lca_count))  %>% # sort descending
  filter(!is.na(last_common)) 
print(prepare)




(lca_plot = ggplot(prepare)  +
    geom_col(aes(x=fct_reorder(last_common, lca_count, .desc=T), y=lca_count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90))
)

# fix names to include only LCA
lca2 = function(x) {
  require(dplyr)
  taxnames = c('superkingdom', 'phylum', 'order', 'family', 'genus', 'species')
  shortnames = apply(x[,taxnames], 2, unique)
  countshnames = sapply(shortnames, length)
  numcount = countshnames==1
  lastuni = tail(names(shortnames[numcount==T]), n=1)
  nombre = as.data.frame(x[1,which(colnames(x) == lastuni)])
  newtax <- as.list(ifelse(countshnames==1,shortnames,NA))
  
  ret = x %>% 
    mutate(last_common = as.character(nombre[[1]])) %>%
    mutate(level_lca = lastuni) %>%
    mutate(superkingdom = newtax$superkingdom) %>%
    mutate(phylum = newtax$phylum) %>%
    mutate(class = newtax$class) %>%
    mutate(order = newtax$order) %>%
    mutate(family = newtax$family) %>%
    mutate(genus = newtax$genus) %>%
    mutate(species = newtax$species)
  return(ret)
}

test = blasthits[500:1000,] %>% 
  group_by(QueryID) %>%
  group_modify(~ lca2(.))

cluster <- new_cluster(8)
cluster_library(cluster, 'dplyr')
cluster_copy(cluster, 'lca2')


#split groups across multiple CPU cores
prepare = blasthits %>%
  group_by(QueryID) %>%
  partition(cluster) %>%  #split groups into cluster units
  do({lca2(.)}) %>%
  collect()

gencount = prepare %>%
  group_by(QueryID) %>%
  slice(1) %>%
  group_by(genus) %>%
  summarize(count=n()) %>%
  arrange(desc(count)) %>%
  filter(!is.na(genus))

(p1 = ggplot(gencount) +
    geom_col(aes(x=fct_reorder(genus, count, .desc=T), y=count)) + 
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab('Genus')
)
