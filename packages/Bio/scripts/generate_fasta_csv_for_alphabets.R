require(tidyverse)
library(random)

alphabetDna <- c('A','C','G','T')
alphabetRna <- c('A','C','G','U')
alphabetPt <- c('G', 'L', 'Y', 'S', 'E', 'Q', 'D', 'N', 'F', 'A',
                'K', 'R', 'H', 'C', 'V', 'P', 'W', 'I', 'M', 'T')

toAlphabet <- function(v, alph){
  paste(sapply(v, function(ci){ alph[ci]; }), collapse = '')
}

seqPt <- toAlphabet(sample.int(20, 35, replace=TRUE), alphabetPt);
seqDna <- toAlphabet(sample.int(4, 35, replace=TRUE), alphabetDna);
seqRna <- toAlphabet(sample.int(4, 35, replace=TRUE), alphabetRna);
# probability to mutate
seq_p1 <- c(100,  100,  100, 100, 100,   5,  2,  2,  50, 3,
           100,  100,   7,  2,  2,   7,  2,  33,  100, 100,
           100,  100,  100, 100, 100, 100, 100, 100, 100, 2,
           100,  100,  100, 100, 100)
seq_p2 <- c(100,  100,   7,  2,  2,   7,  2,  33,  100, 100,
           100,  100,  100, 100, 100, 100, 100, 100, 100, 2,
           100,  100,  100, 100, 100,   5,  2,  2,  50, 3,
           100,  100,  100, 100, 100)

# mutate string s with probability p and alphabet
seq_mutate <- function(s, p, alphabet){
  # s <- seqDna
  # p <- seq_p
  # alphabet <- alphabetDna
  res_s <- s
  res_p <- p
  for (i in 1:(str_length(res_s)*2)) {
    pos <- sample.int(str_length(res_s), 1)
    if (sample.int(100, 1) < res_p[pos]) {
      cast <- sample.int(100, 1) # mutation type probabilty
      if (0 < cast && cast <= 2 ) {
        #insertion
        res_s <- paste(substr(res_s, 1, pos), alphabet[sample.int(4, 1)], substr(res_s, pos+1, str_length(res_s)), collapse='', sep='')
        res_p <- c(res_p[1:pos], c(100), res_p[(pos+1):length(res_p)])
        #cat('insertion');
      } else if (2 < cast && cast <= 4 ) {
        # deletion
        res_s <- paste(substr(res_s, 1, pos-1), substr(res_s, pos+1, str_length(res_s)), collapse = '', sep='')
        res_p <- c(res_p[1: (pos-1)], res_p[(pos+1):length(res_p)])
        #cat('deletion');
      } else {
        # replace
        res_s <- paste(substr(res_s, 1, pos-1), alphabet[sample.int(4, 1)], substr(res_s, pos+1, str_length(res_s)), collapse='', sep='')
        #cat('replace')
      }
      #cat(res, '\n')
    }
  }
  res_s;
}

for (n in c(100,1000,10000, 100000, 1000000)){
  fastaDna_df <- data.frame(id = 1:n, sequence = sapply(1:n, function(id){ seq_mutate(seqDna, seq_p1, alphabetDna)}));
  write_csv(fastaDna_df, sprintf('../files/data/sample_FASTA_DNA-%d.csv', n));
}

fastaRna_df <- data.frame(id = 1:100, sequence = sapply(1:100, function(id){ seq_mutate(seqRna, seq_p2, alphabetRna)}));
write_csv(fastaRna_df, 'D:/HOME/atanas/Datagrok/projs/public/packages/Bio/files/samples/sample_FASTA_RNA.csv');

fastaPt_df <- data.frame(id = 1:100, sequence = sapply(1:100, function(id){ seq_mutate(seqPt, seq_p2, alphabetPt)}));
write_csv(fastaPt_df, 'D:/HOME/atanas/Datagrok/projs/public/packages/Bio/files/samples/sample_FASTA_PT.csv');
