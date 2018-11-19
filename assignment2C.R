# I have taken Mollusca as my biological question and created a phylogenetic hypothesis for it.
#Package installations------
##install.packages("rentrez")
library(rentrez)

##install.packages("seqinr")
library(seqinr)

##install.packages("Biostrings")
library(Biostrings)

##install.packages("stringr")
library(stringr)

##install.packages("tidyverse")
library(tidyverse)

##install.packages("dendextend")
library(dendextend)

##source("https://bioconductor.org/biocLite.R")
##biocLite("Biostrings")
##biocLite("muscle")
##biocLite("DECIPHER")
##biocLite("ape")
library(DECIPHER)
library(muscle)
library(Biostrings)
library(ape)

#to find the list of available databases
entrez_dbs() 
#search for records of the genus Mollusca in the NCBI nucleotide data base.
search <- entrez_search(db = "nuccore", term = "Mollusca[ORGN]", retmax = 100) 

Mollusca_fetch <- entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta")#Download the data in FASTA format and convert it to a data frame.
Mollusca_fetch

# What class is it?
class(Mollusca_fetch)



# Write to file in your current working directory.
# Notice the new line (\n) deliminator.
write(Mollusca_fetch, "Mollusca_fetch.fasta", sep = "\n") 
# Read it back in as DNA StringSet using the readDNAStringSet() function.
stringSet <- readDNAStringSet("Mollusca_fetch.fasta")
#How many total sequences do you have? Provide the line of code and the answer
dfMollusca <- data.frame(Name = names(stringSet), Mollusca_Sequence = paste(stringSet))
dfMollusca$Name <- word(dfMollusca$Name, 2L, 3L)
#I have 100 sequences in total 

#How many unique species do you have? Provide the line of code and the answer.
length(unique(dfMollusca$Name))

dfMollusca_Subset <- dfMollusca %>% 
  group_by(Name) %>%  ## Group by Name.
  sample_n(1)  ## Handy function in dplyr for randomly selecting rows from tbl object.
#I have 3 unique species 

#Perform a multiple sequence alignment (MSA). Provide the code, and explain your choice of alignment algorithm and chosen arguments

Mollusca_alignment1 <- DNAStringSet(muscle::muscle(stringSet, log = "log.tx", verbose = T))

length(Mollusca_alignment1[[1]])
lapply(Mollusca_alignment1, str_count, ("-"))
mean(unlist(lapply(Mollusca_alignment1, str_count, ("-"))))


#Converting data type for further downstream analysis. 
dnaBin_Mollusca <- as.DNAbin(Mollusca_alignment1)
###############IMPORTANT EDIT change names on dendrogram to ids, allows it to be readable, and prevents it from running so long
numbers <- (1:100)
names(dnaBin_Mollusca) <- numbers

#creating a distance matrix to use to create dendrogram
distanceMatrix1 <- dist.dna(dnaBin_Mollusca, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
#Cluster your sequences into OTUs 

clusters_Mollusca <- IdClusters(distanceMatrix1,
                              method = "single",
                              cutoff= 0.02,
                              showPlot = TRUE,
                              type = "both",
                              verbose = TRUE)
#Check: Should be getting a list with a dendrogram and clusters
clusters_Mollusca

#Present a visualization of your clusters 
####The hang on the dendrogram causes a shift in the height of the leaves, the purpose of this change was also clarity.
plot(hang.dendrogram(clusters_Mollusca[[2]], hang = 0.03, ), main = "Cluster dendrogram")
length(unique(unlist(clusters_Mollusca[1][1])))

#Then I have took biological data of Mollusca and created the phylogenetic hypothesis.

search <- entrez_search(db = "nuccore", term = "Antimicrobial peptides derived from molluscs", retmax = 100)
Mollusca_fetch2 <- entrez_fetch(db = "nuccore", id = search$ids, rettype = "fasta")
Mollusca_fetch2

write(Mollusca_fetch2, "Mollusca_fetch2.fasta", sep = "\n") 
# Read it back in as DNA StringSet using the readDNAStringSet() function.
Mollusca_stringSet <- readDNAStringSet("Mollusca_fetch2.fasta")
#How many total sequences do you have? Provide the line of code and the answer
dfMollusca2 <- data.frame(Name = names(Mollusca_stringSet), Mollusca_Sequence = paste(Mollusca_stringSet))

dfMollusca2$Name <- word(dfMollusca2$Name, 2L, 3L)
#I have 100 sequences in total 

#How many unique species do you have? Provide the line of code and the answer.

length(unique(dfMollusca2$Name))

dfMollusca2_Subset <- dfMollusca2 %>% 
  group_by(Name) %>%  ## Group by Name.
  sample_n(1)  ## Handy function in dplyr for randomly selecting rows from tbl object.
#I have 3 unique species 

#Perform a multiple sequence alignment (MSA). Provide the code, and explain your choice of alignment algorithm and chosen arguments

Mollusca_alignment2 <- DNAStringSet(muscle::muscle(Mollusca_stringSet, log = "log.tx", verbose = T))

length(Mollusca_alignment2[[1]])
lapply(Mollusca_alignment2, str_count, ("-"))
mean(unlist(lapply(Mollusca_alignment2, str_count, ("-"))))

#Converting data type for further downstream analysis. 
dnaBin_Mollusca2 <- as.DNAbin(Mollusca_alignment2)

#creating a distance matrix to use to create dendrogram
distanceMatrix2 <- dist.dna(dnaBin_Mollusca2, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
#Cluster your sequences into OTUs 
clusters_Mollusca2 <- IdClusters(distanceMatrix2,
                                method = "single",
                                cutoff= 0.02,
                                showPlot = TRUE,
                                type = "both",
                                verbose = TRUE)
#Check: Should be getting a list with a dendrogram and clusters
clusters_Mollusca2

#Present a visualization of your clusters, again this one has the hang as
plot(hang.dendrogram(clusters_Mollusca2[[2]], hang = 0.03, ), main = "Cluster dendrogram")

###############dend_diff is a function for putting two dendrograms unto the same plot, so it can be compared properly
total <- dend_diff(clusters_Mollusca[[2]], clusters_Mollusca2[[2]])
total


