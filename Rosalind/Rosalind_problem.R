# 1. Count ACGT

dnaseq <- "GCGATCCATGTAAGTACCGCCATGACAGATAGTACACAGGACTGCGTGTAGGTTTAGCTTATGGAACCCAGTCAATTTTTAAAGTGCCAACTCGCGTTGCTTAGCCAGATCTAATTTGAAGGTATAGTTACCTCAGGCCGCGTTGGTTCTTCTCTGCCTTCCTGGGGTTTACAACAAGGCGATGGTCCTGACTTTGGTCGGTACGTTGAGGGGATCAAGCTCTAAGGCTATCGTAATAATTAATGATCGTTCCGCCTTGCTCTGGTTCGTTGTCGATTCATGGTTGTTCGCACTTTGTTCGTCCATGGCGCGCGAACCAAAGAGAATCATCATCCCGCGCTCTCTCAACGTGGGCTACGATACATACCATCAGTGGAAAGCATGGAGCGCATGTGTTCTCGTTTGCTATGTCTAATGAACTGTAGGAGACCTGTGCACTCCTATGCTGGATTGGACATGTACCGACCCCGTCTATTAGGTTTCCGCGGACTCCAATTTAGACCATCGTGCCATGCGACGACTCCACACTTTGGGTAGCCGACTCAGTTTCGGATAACTGTGGAATGAAGTAGAGGGGGGCTGTAGCTACACAAAAGTAAATGCTCCCACACGAATCTCCACGTGCGACGTTCGCTTGTTTATTTTCTGGGCTACCGCATAACGATTGGGGATGCTGCGTTGTTAAAGAAAATCCAGTATGCCGCTGGCGTCTGCAGACGCTTCCTGTGGCACCATATATGCCCCTTTCTTGTGGTAACCGGACAAGGTTCGTGCACCAGTGTCTGATCCCTATTTTTGGAAGACTAGTGATACTGGGCCTTAGGGACTTCCTCGTGGCACCTGTGTCCTGTATCGCGGTGCAAGGAAACGACCCGCCCAAATCGTTACCAAATCCCCAAAATAGGCTGGGTGAAGTCTACTACCTTTTAAAAATTCCTCCACGGTAATACTCGGTTTCGCGTGTATATTCTCCTCCTAAGATGTCTAATTGAACGCCAT"

# string split by base
base <- strsplit(dnaseq, "")
# Check Freq.
table(base[[1]])

# 2. Transcribing DNA into RNA

# DNA to RNA (T -> U)
dnaseq <- "GGGGTCGCCCTGAACTTACAAATGTTACTCAAATCATTGACGGAAATCTTTGAAAGCGCGGACCCAAAGTTAGACCATTTAACCCACCTAAATCGGGAGGGTGTATGCCATTTACTCCAAAGTTACCGCCGGAGTCGGTTGCAAACGGCGATTAATTAATATGGCCCCTATTGGGATAACCTTTTGTAGATACTGCAAGATGGTAGTGACAACAGGATCGGCTAGATGGAACATCCATGGCTACCAGCGGGGAGGAATGGGGCCGCCTGGATCGTCCGGCATTCGTAATACCGACGACGATGCGACTATGTTCATCAGGACCTAACGAGGTGCTACGAAAAGTGACTACGATAAGACGTCTTCGAGCCTGATGCGCGTAGTTTGACTAGCCTGTCGCAGCGCTTGTCAGTGTCCGACATGGGTAAAACATAGCCAACTAGTCGCATGGCATTCGCCATACAAGGCCAGCATAGGTAACAGTGTCAATGTAAGCCGTACCCAGCTGTGGGAGCTGCCCCTGCGCGACAGGTAATGAATTTCCACTCGTTGCTATCCGTAAGTGACCACATTTCATTGAGCCTAGAGCGGATTTTGGTGCGATTTGTCTAACCCCGCATCGCAAACACAAATTTTTAGTGAGGAGCATAATTTGGCACCTTTACAGTATATCCTGTCGCCGACTGGGATAACGAGGTATATAGAAAACAAGCATGCGCGCATTCTACGCTGCGGTTGCTATGCACTAAGTCAGTAGAGACTTTGTATTCCCAGGTTGGCTGGCGACGCGACATGCTCTTCGAAGAGCTCTTCGCAGCCGGTGTCACACTCCCACAATGGCACATTAACTTCTTGCGGCGCAGAAGAGTGTTGGACCTCTGATCCTTTGTCTAGTATGTGACGAATGTACTAAAATGACGCAACTTCATTATGTACGCTCTATTAATACGGACTGGGCTGTTTCACCCCCATTTAAAAGGAGTCCTAA"
gsub("T","U",dnaseq)

# 3. reverse complementary

dnaseq="GATTCATAGTGGTACTGTCCCCTCCACAGTGACTAACGAAGCTTATCAGAGTCATGTGCGGGACGCCTCCTTGGAATCGCTGTTCCCGAATAAAGTCTATTGTAGCCTCCTAGCTGGTTGACCAAAAGAAATAGCTTAGTTAGTTGATGGTGTAGGATTCACAAACCAAAGTTCCCTGCTCACTGTCGTCGCCTGGAGGGGGGTAGACCCGCATGGCATAAACTAGATGTTCGAGGAATGACTCCAGAGCCTGGAAAGGACATACAGATACAAAAGTGTGGTAATATTCGACGGTTTGAGCGCTCCGGATGCCAGTTCGTCTACAGATTTAAATCGGGGCAGTTGTGATGTCAGGTCGTGTTTTGGCCCCATTGCGTTGAGGAATGGACGCTTCGTCTCGAGAGGTGCTACCAACCAGCTGATATCACAAGAGCGGCGACCGCGAGGCCGGTAGCCGGCCTACTGTGACGACTTTCGCCTTGTCCTAAGACGGTTAATCGCAACTTCGAACTATAAGTATCCGCGCTTAGCCCCGCTCTCAACATGTATCCACTCTCTAGATAGTTGGCTGGCCCTAATGGTAGCTTTATGGCACTCCGTAGTGCCGGACCACTGGCGAAAATTGGGGAATGCGTAGGCTTCTGAGCTTAACGTCCTCTGCCCAACACTCTGGATGCGTATCCCATAGCTTGAGACGTCCCAGTTGCTCGCAGACTAGCCTCCGTAGCACTAGGCGGACAGGCTTTATGGGTGAACAGAGGTCGTGGCGTACTCCTATGGGGCCGGATCTGGTATTAACTGTGGGAGGCCCAGTCCAAGTTCAAGGCACTGCATTAAGGGAGGATGGGACGATATAACACCGGACATCCGTGTAGTCACCATCCGGTTGTGCTGTGATTCGCCAA"
# read each bases

bases <- unlist(strsplit(dnaseq,""))
# reverse it
revbases <- rev(bases)
# empty list for complementary seq
n_str <- c()

# reading bases and convert to complementary seq.
for(i in 1:length(revbases)){
  if(revbases[i] == "A"){
    n_str[i] <- "T"
  }
  if(revbases[i] == "T"){
    n_str[i] <- "A"
  }
  if(revbases[i] == "C"){
    n_str[i] <- "G"
  }
  if(revbases[i] == "G"){
    n_str[i] <- "C"
  }
}

# Alternative
################################################################################
xx <- unlist(lapply(revbases, function(x){
  if(x == "A") n_str <- "T"
  if(x == "T") n_str <- "A"
  if(x == "C") n_str <- "G"
  if(x == "G") n_str <- "C"
  return(n_str)}))

#collapse
paste(xx, collapse = "")
################################################################################

# 4. Fibonacci number (reprocuing k rabbit paires, not only 1 pair)

# f1, f2 = 1
# fn = f(n-1) + f(n-2); f(n-2) is equal to new offspring pair in fn (when reproduction is 1; k = 1)
# given number n=<40 (Month), k =<5 (offspring from single pair)

Fibo <- function(n=5, k=3){
  f <- c()
  f[1] <- 1
  f[2] <- 1
  if(n >= 3){
    for (i in 3:n) {
      f[i] = f[i-1] + k*f[i-2]
    }
  }
  return(f[n])

}
Fibo(5,3)
Fibo(1,1)
Fibo(36,2)

# 5. GC content

# ID: start with > and finish \n
# Seq: Start with A|C|G|T and finish
# Data frame IDs, Seqs
# Number of G + C / Total number

library(Biostrings)
fasta <- readDNAStringSet("../rosalind_gc.txt")
fasta_df <- as.data.frame(subseq(fasta))
IDs <- rownames(fasta_df)

# loop for all GC contents

GCcal <- function(IDs){
  GCs <- c()
  for (i in IDs) {
    nG <- unname(table(strsplit(fasta_df[i,],""))["G"])
    nC <- unname(table(strsplit(fasta_df[i,],""))["C"])
    GCs[i] <- (nG+nC)/length(strsplit(fasta_df[i,],"")[[1]])*100
  }
  return(GCs[order(GCs, decreasing = T)][1])
}

GCcal(IDs)

# 6. counting point mutations

# read bases
# count = 0
# if not same count+1
# return count

seq1 = "CCAAGGCGAAATTTCATATCGTGTCCTCTTTTAGGGACGTTTCCCCCCGAGCTAGCGGACCTTAACTTATTCCCGCGAAACTCGCACTCTTTTTTGTGTGTATTTTTGGTAAGAAGGCAAATGGGTTCTACACCGGCGGCCGACGCCGGCTCACCCATATGGCAGGACGTGTAGGGGCGGATGGTATTCTACAAAGGCGACTAACGACGAAAACGGCCAGTGTAGGTAGAGCTACTGCGTCATCGCCCACTTCTAAATGTGCGTTTCCGGCACTCACCCCGGTTCGCGAATGCGGTTCTCTGCACACGGGTGTTAGGTTTATTTTTCAACAATTTCTCAACGTGCTGAAGTAAGAGTCGGTGCCCAGAAGGGCACTCCATGCCCCATTCATCATCACCTTCCGTTATACTAACTTAGGATGCCATGGCAGTTGACGTAATTCGGCGGCTTTCGGGTCTAGTGTTCATGTGATGGTGCACTAGTTATGGGATTATCTTACACATTCTCCTCAGGAGAGACACCACAGCGGGCATTTACCATGAGCGCATCGCACTTCAAGCGAGAACTGCAAATTGATCCAGCCTCTCTTGTAGCCAATGGTCTCCCTCACCGAATTTGTTGGACTTACTACCTTCATGGTATTGCTATGGCCCCCGATGACGATAACGATGCTAACAATAATGGCCGCTGATGTACGTTAAATACAATCTACTACGTTATCATTTATATCCAGCTAACATTCGGAAATGCTTGTAGCGCTATTAGTATCCATCGCCTTTCGGGCCTACTCACGTCCAACAGTGAATTTGTCAGTGGAGGTTGGGCATTAAGTGGTCAGGTTTATACGTGCAGTTATGTTGTCATGATGTTTTCGTCGGGTCTAATCAAACTTTTTAACGAGAATAGATGCAGATTCGTTGATTAGCACATTCCACTATGGACTATAGACCCACTGAGACCCTGTAATAGGCTTTGCAGCTTTGAAGCCT"
seq2 = "CTAGAGAGGAATTGATAAAGGGTTATGGTTAAAGCGTCTATTGTCAATCAGCGCGCTCTTCACGTGATGTCTGTCCGAAAGTCCCTTTCGAATTTATGCTCATTGATTGCAACAGAGCCAGTGGCATGTACAACGGTTGTCTACCCAGAACTAACTCGATGGCCGATCGTGTAGGTACTAATGCCATTTAAGTCAGGATACCTAGGAGATTATGATTAAGAGTCCGTAGAGCTAGAGCATACTCAGTCCACTCTAGGGGTCTCCTGCTTACTCGCGGTCTGGGCGACGGTTGGGCTTTTCTTTATGTGCGGTTTATGGTCCATTTTCATATTTGCCACAAGAGGGTGACGTTTCGCGCCTTGTTCAAAAGGGAACGCGAGGGTCCTATTATCTTCACGTTCCGGTGGTGTATGAGAAGTTTCCCGTGGTCGGGAAATGTTGCGTAGCTTTTTGAGCTTACAGTACAAGCGAGGATGAACTGCGGCTGTCGATGAGTCAGACACTTCCGTCTAACTCGAGTTCTGATTGAGCATTAATCAAGAGCGCATACCATGTCGGTAAAATTTTGCACATACATCCTACCAACATAGTTGTCAATCGTTCGAAACACGATTTTTGTTTGACCGCCTGCCTTCAGGGTATTGGTGGGCCTGGCACTTGTATTTACGCTGCACACCATATTTGCCACTCTAGTCAGTGGACAATCAGATAATACATGCTGATATGAAACCAACGATTATGAGCCAACGCAGGAAGCACTATACGTATCTATCGAGGACCGGCGTGACAAAAATCCAAGCAAAAATATTCGCAAAGAGAAAAGATATGAAGCCGAGAGATTTAATTGGGCGGTTAAGCTACGGTGTGATTTGCGTAGGGTTCAAGAAAAGTATGGTCCGATGATTGCGGGCGATTCATGATTCACCACAACCAACTTGTGGGTATCTTACACCTGTGAGAGGTTCATATCCTTAAAAAGTTTTAAGGTT"

seq1 = strsplit(seq1,"")[[1]]
seq2 = strsplit(seq2,"")[[1]]

pmc <- function(seq1, seq2){
  counts <- 0
  for (i in 1:length(seq1)) {
    if(!seq1[i]==seq2[i]){
      counts <- counts+1
    }
  }
  return(counts)
}

pmc(seq1, seq2)

# 7. Mendel's First Law

# Probability calculation
# 3 different organisms k, m ,n
# k is dominant Homozygous XX
# m is heterozygous Xx
# n is recessive homozygous xx
# When mating, calculate probability that offspring is dominant phenotype (XX or Xx)
# 1-Pr(recessive)
# Given numbers are number of k, m, n and total population is k+m+n

# Pr(recessive)
#M-M form M-M recessive is 1/4
#M-N (or N-M) from M-N recessive is 1/2
#N-N from N-N recessive is 1

mendel <- function(k=1,m=1,n=1){
  total = sum(k,m,n)
  # Pr(recessive mm)
  mm <- m/total*((m-1)/(total-1))*1/4
  # Pr(recessive mn)
  mn <- m/total*n/(total-1)*1/2+n/total*m/(total-1)*1/2
  # Pr(recessive nn)
  nn <- n/total*(n-1)/(total-1)

  Pr <- 1-(mm+mn+nn)
  return(Pr)
}

mendel(25, 27, 21)

# 8. RNA to protein

# need codon data.frame rownames = codon, col = aa
# start i = 1 , 4, 7... end of seq, read i:i+2 convert to aa.
# when find stop codon, stop to read.

rnaseq <- readChar("../rosalind_prot.txt", file.info("../rosalind_prot.txt")$size)
rnaseq <- gsub("\n","",rnaseq)
# Setting up #
# 1. rnaseq each bases convert to list
rnaseq <- unlist(strsplit(rnaseq,"")) # each bases convert to list
# 2. import codon table as data.frame
codon.df <- read.delim("../codon_table copy.csv", header = T, row.names = 1)
codon.v <- codon.df$Letter
# 3. generate codon named vector
names(codon.v) <- rownames(codon.df)

i <- 1 # base position
protein <- c() # protein seq

# While position of base is smaller than -2 end posotion
while (i < length(rnaseq)-2) {
  codon <- paste(rnaseq[i:(i+2)], collapse = "") # read each 3 base -> codon
  if(!codon.v[codon] == "O"){ #If it is not stop codon
    protein <- c(protein,codon.v[codon])
    i <- i+3
  }
}
protein <- paste(protein, collapse = "")

protein

# 9. find motif
# find pattern and location
# while i< length of seq - length of pattern
# read bases seq[i:(i+length of pattern)]
# if seq[] == motif : record i

seq = "AATCCTACATCCTACATATCCTACGCATCCTACCGGATCCTACAGATCCTACATCCTACAACGATCCTACCTCCTAATCCTACATCCTACGATCCTACCGATCCTACATGCTTCGCAATCCTACAATCCTACATCCTACCCATCCTACATCCTACTCGAATCCTACTTATCCTACATCCTACAAGGTATCCTACCGTTAGGACGAAAATCCTACTGGCATCCTACTATCCTACGCGTATCCTACCGCGCGATCCTACATCCTACGCCGTCATCCTACTAATCCTACCAATCCTACATCCTACCATCCTACCCTGTATCCTACATCCACGTAATCCTACTACCAACTATCCTACATCCTACTAAGCATCCTACTCGACTGCTTGAGATCCTACAAATCCTACAATCCTACTATCCTACATCCTACATCCTACTTGATCCTACATCCTACATCCTACATCCTACATCCTACCGATCCTACCATCCTACATCCTACGTAAGATCCTACGACGCATCCTACATCCTACTATCCTACCGATCCTACCTATCCTACGATCCTACATCCTACAAATCAGACGATCCTACCACCGCATCCTACAATCCTACATATCCTACTGATCCTACATGGCATCCTACAATCCTACCGCGGTATCCTACCGCCGGGCAAATATCCTACGCGCAATCCTACATCCTACATCCTACCATCCTACAACGAATCCTACGTATCCTACATCCTACATCCTACGGAAATCCTACCATCCTACCATCCTACATATCCTACTAAACGAGCTTATCCTACATCCTACGATCCTACTGCTATCCTACACCGGTCCGTGTGAT"
motif = "ATCCTACAT"

seq = unlist(strsplit(seq,""))
location <- c()
i <- 1
while(i<(length(seq)-nchar(motif)+2)){
  subseq <- paste(seq[i:(i+nchar(motif)-1)], collapse = "")
  if(subseq == motif){
    location <- c(location, i)
  }
  i <- i+1
}

location

# 10. Consensus and Profile

# given seq matrix (m x n)
# Freq of each column
# read fasta files by Biostrings::readDNAStringSet()
library(Biostrings)
seq.table <- readDNAStringSet("../rosalind_cons.txt")
# take only seq dataframe
seq.df <- as.data.frame(subseq(seq.table))
# generate matrix by stringr::str_split_fixed(df$selected_col, "", length of seq)
library(stringr)
seq.mx <- str_split_fixed(seq.df$x,"",nchar(seq.df$x[1]))

profile <- data.frame(row.names = c("A","C","G","T"))
i <- 1
consensus <- c()
# i = 1:nchar(seq.df$x[1])
for (i in 1:nchar(seq.df$x[1])) {
  # read freq of each bases A,C,G,T and put in to profile data frame
  profile[1,i] <- table(seq.mx[,i])["A"]
  profile[2,i] <- table(seq.mx[,i])["C"]
  profile[3,i] <- table(seq.mx[,i])["G"]
  profile[4,i] <- table(seq.mx[,i])["T"]
  # named vector for each column
  cons <- c(profile[,i])
  names(cons) <- c("A","C","G","T")
  # find max value to check base, put in cosensus list and collapse
  consensus <- paste(c(consensus,names(cons[order(cons, decreasing = T)[1]])), collapse = "")
}
# NA to 0
profile[is.na(profile)] <- 0

consensus
profile

# 11. Mortal Fibonacci Rabbits

# n-th month if all rabbits live for m months.
# before m month F(n) = F(n-1) + F(n-2)
# on m+1 and m+2 month F(n) = F(n-1) + F(n-2) - 1 (rabbits from first 2 month are eliminated here)
# after m+2 month F(n) = F(n-1) + F(n-2) - F(n- m -1) (rabbits born m-1 month ago are eliminated)
# in R it has a problem because of large integer. But code is correct
# library(gmp)

Fibo <- function(m=80, n=19){
  f=c()
  f[1]=1
  f[2]=1
  i = 3
  while (i <= m) {
    if(i<=n){
      #f[i] = as.numeric(as.bigz(add.bigz(f[i-1], f[i-2])))
      f[i] = f[i-1]+f[i-2]
    } else if(i==n+1 || i==n+2){
      #f[i] = as.numeric(as.bigz(add.bigz(f[i-1],f[i-2]) - as.bigz(1)))
      f[i] = f[i-1]+f[i-2]-1
    } else {
      #f[i] = as.numeric(as.bigz(add.bigz(f[i-1],f[i-2]))-as.bigz(f[i-n-1]))
      f[i] = f[i-1]+f[i-2]-f[i-n-1]
    }
    i=i+1
  }
  return(f[m])
}
Fibo(80,19)
format(Fibo(80,19), scientific = F)

# 12. Overlap Graphs

# k is length of pre- and suffix for DNA string
# Given seq and k = 3, The adjacency list corresponding to k=3. You may return edges in any order

# read fasta file -> data.frame
# convert to list
# k length suffix string[n]== k length prefix string[m], then shows name of string[n] and name of string[m]

library(Biostrings)
seq.table <- readDNAStringSet("../rosalind_grph.txt")
seq.df <- as.data.frame(subseq(seq.table))

overlapgraph <- function(seq.df=seq.df, k=3){
  suffix <- c()
  prefix <- c()
  # add suffix and prefix to seq.df
  for (i in 1:nrow(seq.df)) {
    suffix <- c(suffix,paste(unlist(strsplit(seq.df$x[i],""))[(nchar(seq.df$x[i])-k+1):nchar(seq.df$x[i])], collapse = ""))
    prefix <- c(prefix,paste(unlist(strsplit(seq.df$x[i],""))[1:k], collapse = ""))
  }
  seq.df$suffix <- suffix
  seq.df$prefix <- prefix

  results <- ""
  for (i in 1:nrow(seq.df)) {
    # take samples: suffix of seq[i] == prefix of seq[all]
    nm <- rownames(seq.df)[which(seq.df$suffix[i] == seq.df$prefix)]
    # remove self looping sample
    nm <- nm[!nm %in% rownames(seq.df)[i]]
    # if samples are not empty
    if(!length(nm)==0){
      for (j in nm) {
        # add to resutls suffix string[n] == prefix string[m]
        results <- paste(results,rownames(seq.df)[i], j,"\n")
      }
    }
  }
  # write table without rowname colname quote
  write.table(results, "../results.txt", row.names = F, col.names = F,quote =F)
}
overlapgraph(seq.df, 3)

# 13. Calculating Expected Offspring
# Probability of Dominant offspring
# first Pr(X), second round Pr(X) E(X+X) = Pr(X)*2
# AA - AA : Pr(dominant) is 1 any of allele could result in dominant
# AA - Aa : Pr(dominant) is 1
# AA - aa : Pr(dominant) is 1
# Aa - Aa : Pr(dominant) is 0.75
# Aa - aa : Pr(dominant) is 0.5
# aa - aa : Pr(dominant) is 0

dom_off <- function(input="17326 16933 17814 18799 18625 16430"){
  num <- strsplit(input, " ")[[1]]
  dom <- (as.numeric(num[1])+as.numeric(num[2])+as.numeric(num[3])+as.numeric(num[4])*0.75+as.numeric(num[5])*0.5)*2
  return(dom)
}
dom_off()

# 14. Finding a Shared Motif
# Note. better to use binary search algorithm.
# DNA string 1, length = n
# possible length of substring 2:n
# read from beginning (start from length 2 to n) compare all other seq
# exchange longest common motif

library(Biostrings)
seq.table <- readDNAStringSet("../rosalind_lcsm.txt")
seq.df <- as.data.frame(subseq(seq.table))
n=nchar(seq.df$x[1])
comm_motif=""

# k = 1:(nchar(seq.df$x[1])-1) k is minimum length of motif-1
# i is position 1:(n-k)
for (k in 1:(n-1)) {
  for (i in 1:(n-k)) {
    chars <- paste(unlist(strsplit(seq.df$x[1],""))[i:(i+k)], collapse = "")
    if(!"FALSE" %in% as.character(grepl(chars, seq.df$x, fixed = TRUE)[-1])){
      comm_motif <- chars
    }
  }
}
comm_motif

# it takes more than 5 min, so I need different way
# below from Li Yutze
fn="../rosalind_lcsm.txt"
sharedSubstringUp <- function(fn) {
  g <- as.character(ape::read.FASTA(fn))
  g <- toupper(sapply(g, paste, collapse = ''))
  base <- c('A','C','G','T')
  candidates <- base
  repeat {
    ## remove candidates that are not presented in all sequences
    present <- sapply(candidates, function(subseq) {
      all(grepl(subseq, g))
    })
    if (!any(present)) break
    candidates <- candidates[present]
    ## append 4 base to each remained candidates
    enlong <- unlist(lapply(candidates, paste, base, sep = ''))
    present <- sapply(enlong, function(subseq) {
      subseq <- sub('^.', '', subseq)
      any(subseq == candidates)
    })
    last <- candidates
    candidates <- enlong[present]
  }
  cat(last, '\n')
}
sharedSubstringUp("../rosalind_lcsm.txt")

# binary search algorithm

sharedSubstringDown <- function(fn) {
  g <- as.character(ape::read.FASTA(fn))
  g <- toupper(sapply(g, paste, collapse = ''))
  lengths <- sapply(g, stringr::str_length)
  shortest <- strsplit(g[which(lengths == min(lengths))[1]], '')[[1]]
  g <- g[which(lengths == min(lengths))[-1]]
  len_shortest <- length(shortest)

  l <- 1
  r <- len_shortest
  while (l + 1 < r) {
    cut <- floor({l + r} / 2)
    any_match <- F
    for (i in 1:{len_shortest - cut + 1}) {
      sub <- paste(shortest[i:(i + cut - 1)], collapse = '')
      if (all(grepl(sub, g))) {
        any_match <- T
        break
      }
    }
    if (any_match) l <- cut
    else r <- cut
  }
  cat(sub, '\n', sep = '')
}

# 15.

1*1%%1000000 # possible M
4*1%%1000000 # possible A
3 #possible stop codon

# 16.
IDs <- "Q67JS9
P19246_NFH_MOUSE
P0A4Y7
P01047_KNL2_BOVIN
P03415_VME1_CVMA5
Q13VE3
A4TEW1
P02725_GLP_PIG
P08514_ITAB_HUMAN
P01215_GLHA_HUMAN
P80069_A45K_MYCBO
Q6A9W5
Q32LI2
B4U0J5"

protein_motif <- function(IDs){
  IDs <- unlist(strsplit(IDs, "\n")[[1]])
  motif <- c("N","[^P]","[ST]","[^P]")
  for (ID in IDs) {
    seq.table <- as.character(read.FASTA(paste0("http://www.uniprot.org/uniprot/",ID,".fasta")))
    seq.df <- toupper(sapply(seq.table, paste, collapse = ""))
    seq = unlist(strsplit(seq.df,""))
    location <- c()
    i <- 1
    while(i<(length(seq)-length(motif)+2)){
      subseq <- paste(seq[i:(i+length(motif)-1)], collapse = "")
      if(grepl(paste(motif, collapse = ""),subseq)){
        location <- c(location, i)
      }
      i <- i+1
    }
    cat(paste(ID,"\n"))
    cat(paste(location))
  }
}


# Needleman-wunsch alignment algorithm


needles <- function(pattern, subject, params=defaultNeedleParams) {
  MATCH <- params$MATCH
  MISMATCH <- params$MISMATCH
  GAP <- params$GAP
  GAPCHAR <- params$GAPCHAR

  patt <- strsplit(pattern, "")[[1]]
  subj <- strsplit(subject, "")[[1]]
  ## Initialize
  scoreMatrix <- matrix(NA,     ncol=1+length(patt), nrow=1+length(subj))
  direcMatrix <- matrix("none", ncol=1+length(patt), nrow=1+length(subj))
  scoreMatrix[1,1] <- 0
  for (j in 1:length(patt)) {
    scoreMatrix[1, j+1] <- GAP*j
    direcMatrix[1, j+1] <- "left"
  }
  for (i in 1:length(subj)) {
    scoreMatrix[i+1,1] <- GAP*i
    direcMatrix[i+1,1] <- "up"
  }

  ## Fill
  for (i in 1:length(subj)) {
    for (j in 1:length(patt)) {
      ## Translating from 0-based arrays and vectors to 1-based
      I <- i + 1
      J <- j + 1
      ## Calculate (mis)match scores
      letter1 <- patt[J-1]
      letter2 <- subj[I-1]
      if(letter1 == letter2) {
        diagonalScore <- scoreMatrix[I-1, J-1] + MATCH
      } else {
        diagonalScore <- scoreMatrix[I-1, J-1] + MISMATCH
      }
      ## Calculate gap scores
      upScore   <- scoreMatrix[I-1, J] + GAP
      leftScore <- scoreMatrix[I, J-1] + GAP
      ## Choose best score
      if (diagonalScore >= upScore) {
        if (diagonalScore >= leftScore) {
          scoreMatrix[I, J] <- diagonalScore
          direcMatrix[I, J] <- "diagonal";
        } else {
          scoreMatrix[I, J] <- leftScore
          direcMatrix[I, J] <- "left";
        }
      } else {
        if (upScore >= leftScore) {
          scoreMatrix[I, J] <- upScore
          direcMatrix[I, J] <- "up";
        } else {
          scoreMatrix[I, J] <- leftScore
          direcMatrix[I, J] <- "left";
        }
      }
    }
  }
  theScore <- scoreMatrix[I, J]

  ## backtrace
  J <- length(patt)+1
  I <- length(subj)+1
  align1 <- align2 <- c()
  while(1) {
    direc <- direcMatrix[I, J]
    if (direc == 'none') {
      break
    }
    if (direc == 'diagonal') {
      align1 <- c(patt[J-1], align1)
      align2 <- c(subj[I-1], align2)
      I <- I-1
      J <- J-1
    } else if (direc == 'left') {
      align1 <- c(patt[J-1], align1)
      align2 <- c(GAPCHAR, align2)
      J <- J-1
    } else if (direc == 'up') {
      align1 <- c(GAPCHAR, align1)
      align2 <- c(subj[I-1], align2)
      I <- I-1
    } else {
      stop("This is not supposed to happen.")
    }
  }
  list(score=theScore,
       align1=paste(align1, collapse=''),
       align2=paste(align2, collapse=''),
       sm=scoreMatrix,
       dm=direcMatrix)
}


##-----------------------------------------------------------------------------
needleScores <- function(pattern, subjects, params=defaultNeedleParams) {
  scores <- sapply(subjects,
                   function(x, y, p) {
                     needles(x, y, p)$score
                   },
                   y=pattern,
                   p=params)
  scores
}



accession_ids <- "JX205496.1 JX469991.1"
accession_ids <- unlist(strsplit(accession_ids, " "))
xx <- ape::read.GenBank(accession_ids)
ape::write.dna(xx,paste0("../needleman.fasta"),format = "fasta")
seq1 <- as.character(ape::read.FASTA("../needleman.fasta"))[[1]]
seq2 <- as.character(ape::read.FASTA("../needleman.fasta"))[[2]]
seq1 <- toupper(paste(seq1, collapse = ""))
seq2 <- toupper(paste(seq2, collapse = ""))

library(Biostrings)
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
pairwiseAlignment(seq1, seq1, type="global", gapOpening = 10, gapExtension = 1)

# Transitions and Transversions
# Transition is purine to purine A <-> G || pyrimidine to pyrimidine C <-> T
# Transversion is purine to pyrimidine A <-> C || A<->T || G<->C || G<->T

# given seq1, seq2 transition/transverion Ratio(s1, s2)
# compare each base
# If seq1[i] == seq[i], i+1
# When seq1[i] == A or G, seq2[i] == G or A ; transition +1
# When seq1[i] == C or T, seq2[i] == T or C ; transversion +1

# ratio = transition/transversion

# read Fasta
path = "../rosalind_tran.txt"
ts_tv_ratio <- function(path){
  seq <- as.character(ape::read.FASTA(path))
  seq1 <- toupper(seq[[1]])
  seq2 <- toupper(seq[[2]])
  stopifnot(length(seq1)==length(seq2))
  transition <- 0
  transversion <- 0
  for (i in 1:length(seq1)) {
    if(seq1[i] == seq2[i]){
      next
    } else if(all(grepl(seq1[i],"[A]"), grepl(seq2[i],"[G]")) || all(grepl(seq1[i],"[C]"), grepl(seq2[i],"[T]")) || all(grepl(seq1[i],"[G]"), grepl(seq2[i],"[A]")) || all(grepl(seq1[i],"[T]"), grepl(seq2[i],"[C]"))) {
      transition <- transition+1
    } else {
      transversion <- transversion +1
    }
  }
  return(transition/transversion)
}

ts_tv_ratio(path)


# RNA Splicing given fasta file;

# 1st str is RNA string (convert any of T to U; gsub("T", "U", RNAstring))
# 2nd ... nth are introns
# delete all intron gsub(intron, "", RNAstirng)
# translate to protin seq
path <- "../rosalind_splc.txt"
RNAsplicing <- function(path){
  # read fasta and seperate RNAseq and intron
  # RNAseq
  RNAseq <- as.character(ape::read.FASTA(path))
  RNAstring <- toupper(unlist(RNAseq[[1]]))
  RNAstring <- gsub("T","U",RNAstring)
  RNAstring <- paste(RNAstring, collapse = "")

  # intron
  introns <- lapply(RNAseq[2:length(RNAseq)], function(x){
    paste(gsub("T","U",toupper(unlist(x, use.names = F))), collapse = "")
  })
  introns <- unlist(introns)
  introns <- paste(introns, collapse = "|")

  # substrate introns from RNAseq = exon
  RNAstring <- gsub(introns, "", RNAstring)

  # open codon table
  codon.df <- read.delim("../codon_table copy.csv", header = T, row.names = 1)
  codon.v <- codon.df$Letter
  # generate codon named vector
  names(codon.v) <- rownames(codon.df)

  # preparation of loop
  i <- 1 # base position
  protein <- c() # protein seq

  # While position of base is smaller than -2 end posotion because last 2 position is not start position of codon
  while (i < nchar(RNAstring)-2) {
    rnaseq <- unlist(strsplit(RNAstring, "")[[1]])
    codon <- paste(rnaseq[i:(i+2)], collapse = "") # read each 3 base -> codon
    if(!codon.v[codon] == "O"){ #If it is not stop codon
      protein <- c(protein,codon.v[codon])
      i <- i+3
    }
  }
  protein <- paste(protein, collapse = "")
  return(protein)
}

RNAsplicing(path)

#Calculating Protein Mass


protein_mass <- function(path){
  # read protein seq file
  seq <- readr::read_file(path)
  seq <- gsub("\n","", seq)
  # read and import protein monoisotopic mass as named vector
  mass_table <- read.table("../monoisotopic_mass_protein.csv", sep = " ", header = F, row.names = 1)
  mass_v <- as.numeric(mass_table$V2)
  names(mass_v) <- rownames(mass_table)
  # initializing
  mass <- 0
  seq_l <- unlist(strsplit(seq, "")[[1]])
  # calculate protein mass
  for (i in 1:nchar(seq)) {
    mass <- mass+mass_v[seq_l[i]]
  }
  return(mass)
}

path="../rosalind_prtm.txt"
protein_mass(path)
