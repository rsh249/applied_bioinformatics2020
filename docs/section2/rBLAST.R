library(Biostrings)
library(taxonomizr)
library(ggplot2)
library(ShortRead)
library(dplyr)
library(stringi)
library(forcats)
library(rBLAST)


#devtools::install_github("mhahsler/rBLAST")

setwd('/projectnb/ct-shbioinf')
#fix RStudio PATH issues for blast+ and sratoolkit
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/blast+/2.7.1/install/bin", sep=":"))
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/sratoolkit/2.9.2/install/bin/", sep=":"))

# Download SRA File:
srr=c('SRR11043475')
system(paste('fastq-dump', srr, sep=' '))


# Read taxonomy database
taxaNodes<-read.nodes.sql("./taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("./taxonomy/data/names.dmp")


# read fastq
dna = readFastq('.', pattern=srr)
reads = sread(dna, id=id(dna)) 
qscores = quality(dna) 

# plot readlength
widths = as.data.frame(reads@ranges@width)
(widthplot <- ggplot(widths) +
    geom_histogram(aes(x=reads@ranges@width), binwidth = 10) + 
    theme_linedraw() + 
    xlab('Read Length (bp)') +
    xlim(0,2000) +
    ggtitle('Read length distribution for 450bp amplicon'))

# plot qscores
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across 
avgscores = as.data.frame(avgscores)
(ggplot(avgscores) +
    geom_histogram(aes(x=avgscores), binwidth=0.2) +
    theme_linedraw() +
    xlab('Quality Score') +
    ggtitle('Per Read Average Quality'))

#cluster reads
#stringDist(reads)
#Clust <- hclust(stringDist(reads), method = "single")
#Groups <- cutree(Clust, h = 2)
#head(sort(table(lane2.1Groups), decreasing = TRUE))

# remove duplicates: library(stringi)
u.reads = stri_unique(reads) # Often does not remove anything! NO exact duplicates
# Filter short reads:
#long.u.reads = u.reads[which(widths>400)]
#DNAset = DNAStringSet(long.u.reads)
DNAset=DNAStringSet(u.reads)
#dist=stringDist(reads)
# Try clustering reads by similarity.


## blast
bl <- blast(db="blast/nt.fa")
cl <- predict(bl, DNAset, BLAST_args = '-num_threads 12 -evalue 1e-100')
accid = as.character(cl$SubjectID) # accession IDs of BLAST hits

# Plot results

#takes accession number and gets the taxonomic ID
ids<-accessionToTaxa(accid, './taxonomy/data/accessionTaxa.sql')
#taxlist displays the taxonomic names from each ID #
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)

cltax=cbind(cl,taxlist)

cltop = cltax %>% 
  group_by(QueryID) %>% 
  top_n(10, Bits) %>%
  filter(Bits>750) %>%
  filter(Perc.Ident>85)

ggplot(cltop) + geom_density(aes(x=Alignment.Length))

(ggplot(data=cltop) +
    geom_bar(aes(x=fct_infreq(genus))) +
    theme_minimal() +
    theme(    
      axis.text.x  = element_text(angle = 45, hjust=1)
      ) +
    xlab('')
  )



(ggplot(data=cltop) +
    geom_density2d(aes(x=Bits, y=Alignment.Length)))
