library(rentrez)

by = 1000
genesearch <- entrez_search(db="nuccore", term="rbcL[Gene] AND plants[filter]", use_history=TRUE)

get_gene = function(x) {
  recs <- entrez_fetch(db="nuccore", web_history=genesearch$web_history,
                       rettype="fasta", retmax=by, retstart=x)
  cat(recs, file="rbcL.fasta", append=TRUE)
  cat(x+by, "sequences downloaded\r")
}

sapp = lapply(seq(1, genesearch$count, by), get_gene)
