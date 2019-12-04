setwd("/mnt/Files_1/sprosser/R_Working_Directory") # change working directory to whatever you want

#load applicable libraries (install these packages upon first use)
library(Biostrings)
library(muscle)
library(DECIPHER)

#NOTE: This script is intedced to be run as part of a larger bash script. It will require some minor modifications if run independently. E.g. This script removes
#the string _TARGETS.fasta from each FASTA file name, because the bash script generates files with that suffix. If FASTA files do not have that name, the
#code must be modified accordingly.

#This script will open each FASTA file in working directory and map the reads to their correct position (by identifyin primer sequences in the reads), then generate a 
#consensus sequence for each FASTA file and combine them into a single output FASTA file. It also outputs some statistics as csv files.

#User must provide FASTA files as well as a tab separated text file called "samples.txt" containing the following headers:
#Well ProcessID Order Identification

#Well - This can be anything. We use plate records (A01-H12) to generate this file so well locators make it easy to visualize.
#ProcessID - This is the sample ID.
#Order - Order-level taxonomy of the sample.
#Identification - Lowest taxonomy ID (if genus-species, use underscore between genus and species name)

#NOTE: samples.txt file cannot contain spaces. Replace with underscores if necessary.

#Primer sequences from 5' -> 3'. Change these if other primers are used.
PrimerF1 <- "ATTCAACCAATCATAAAGATATTGG"
PrimerF2 <- "attrrwratgatcaartwtataat"
PrimerF3 <- "ttataattggdggrtttggwaattg"
PrimerF4 <- "TAAGATTTTGANTNYTNCCNCC"
PrimerF5 <- "atttttwswctwcatwtdgcwgg"
PrimerF6 <- "tatttgtwtgakcwrtwkkwattac"
PrimerR1 <- "WGGTATWACTATRAARAAAATTAT"
PrimerR2 <- "TCARAAWCTWATRTTRTTTADWCG"
PrimerR3 <- "ARDGGDGGRTAWACWGTTCAWCC"
PrimerR4 <- "GTWGWAATRAARTTDATWGCWCC"
PrimerR5 <- "AGCTCCWGCTARNACNGG"
PrimerR6 <- "TAAACTTCTGGATGTCCAAAAAATCA"

PrimerR1.rc <- as.character(reverseComplement(DNAStringSet(PrimerR1)))
PrimerR2.rc <- as.character(reverseComplement(DNAStringSet(PrimerR2)))
PrimerR3.rc <- as.character(reverseComplement(DNAStringSet(PrimerR3)))
PrimerR4.rc <- as.character(reverseComplement(DNAStringSet(PrimerR4)))
PrimerR5.rc <- as.character(reverseComplement(DNAStringSet(PrimerR5)))
PrimerR6.rc <- as.character(reverseComplement(DNAStringSet(PrimerR6)))


#This calculates trim legnths based on PacBio UMIs (16 bp) and adapters (30 bp) plus primer lengths.
F1.5p <- 16 + 30 + width(PrimerF1)
F1.3p <- 16 + 30 + width(PrimerR1)
F2.5p <- 16 + 30 + width(PrimerF2)
F2.3p <- 16 + 30 + width(PrimerR2)
F3.5p <- 16 + 30 + width(PrimerF3)
F3.3p <- 16 + 30 + width(PrimerR3)
F4.5p <- 16 + 30 + width(PrimerF4)
F4.3p <- 16 + 30 + width(PrimerR4)
F5.5p <- 16 + 30 + width(PrimerF5)
F5.3p <- 16 + 30 + width(PrimerR5)
F6.5p <- 16 + 30 + width(PrimerF6)
F6.3p <- 16 + 30 + width(PrimerR6)

#This calcualtes the search window in which to look for primer sequences.
Search.Distance <- max(c(F1.5p,F1.3p,F2.5p,F2.3p,F3.5p,F3.3p,F4.5p,F4.3p,F5.5p,F5.3p,F6.5p,F6.3p))+10

#initialize some variables
output <- DNAStringSet()
Reads.Table <- data.frame()
Duplex.Table <- data.frame()

#Read in samples.txt file
Sample.Taxonomy.Table <- read.table(file = "samples.txt",header = TRUE,sep = "\t")

#the following loop will execute for each FASTA file
for (i in list.files(pattern = "*.fasta")){
  #initialize file-specific variables
  Input <- readDNAStringSet(i)
  
  #Extract sample ID from FASTA file name
  if (any(grep(pattern = "_TARGETS.fasta", i))){
    Sample.Name <- gsub("_TARGETS.fasta","",i)
  }
  if (any(grep(pattern = "_READS.fasta", i))){
    Sample.Name <- gsub("_READS.fasta","",i)
  }
  
  #Initialize more file-specific variables
  F1.Holding <- DNAStringSet()
  F2.Holding <- DNAStringSet()
  F3.Holding <- DNAStringSet()
  F4.Holding <- DNAStringSet()
  F5.Holding <- DNAStringSet()
  F6.Holding <- DNAStringSet()
  R1.Holding <- DNAStringSet()
  R2.Holding <- DNAStringSet()
  R3.Holding <- DNAStringSet()
  R4.Holding <- DNAStringSet()
  R5.Holding <- DNAStringSet()
  R6.Holding <- DNAStringSet()
  
  F1.Leftovers <- DNAStringSet()
  F2.Leftovers <- DNAStringSet()
  F3.Leftovers <- DNAStringSet()
  F4.Leftovers <- DNAStringSet()
  F5.Leftovers <- DNAStringSet()
  F6.Leftovers <- DNAStringSet()
  R1.Leftovers <- DNAStringSet()
  R2.Leftovers <- DNAStringSet()
  R3.Leftovers <- DNAStringSet()
  R4.Leftovers <- DNAStringSet()
  R5.Leftovers <- DNAStringSet()
  R6.Leftovers <- DNAStringSet()
  
  F1.Reads.Single <- DNAStringSet()
  F2.Reads.Single <- DNAStringSet()
  F3.Reads.Single <- DNAStringSet()
  F4.Reads.Single <- DNAStringSet()
  F5.Reads.Single <- DNAStringSet()
  F6.Reads.Single <- DNAStringSet()
  
  F1.Reads.Duplex <- DNAStringSet()
  F2.Reads.Duplex <- DNAStringSet()
  F3.Reads.Duplex <- DNAStringSet()
  F4.Reads.Duplex <- DNAStringSet()
  F5.Reads.Duplex <- DNAStringSet()
  F6.Reads.Duplex <- DNAStringSet()
  
  F1.Reads.Trimmed.Single <- DNAStringSet()
  F2.Reads.Trimmed.Single <- DNAStringSet()
  F3.Reads.Trimmed.Single <- DNAStringSet()
  F4.Reads.Trimmed.Single <- DNAStringSet()
  F5.Reads.Trimmed.Single <- DNAStringSet()
  F6.Reads.Trimmed.Single <- DNAStringSet()
  
  F1.Reads.Trimmed.Duplex <- DNAStringSet()
  F2.Reads.Trimmed.Duplex <- DNAStringSet()
  F3.Reads.Trimmed.Duplex <- DNAStringSet()
  F4.Reads.Trimmed.Duplex <- DNAStringSet()
  F5.Reads.Trimmed.Duplex <- DNAStringSet()
  F6.Reads.Trimmed.Duplex <- DNAStringSet()
  
  F1.Reads.Aligned.Single <- DNAStringSet()
  F2.Reads.Aligned.Single <- DNAStringSet()
  F3.Reads.Aligned.Single <- DNAStringSet()
  F4.Reads.Aligned.Single <- DNAStringSet()
  F5.Reads.Aligned.Single <- DNAStringSet()
  F6.Reads.Aligned.Single <- DNAStringSet()
  
  F1.Reads.Aligned.Duplex <- DNAStringSet()
  F2.Reads.Aligned.Duplex <- DNAStringSet()
  F3.Reads.Aligned.Duplex <- DNAStringSet()
  F4.Reads.Aligned.Duplex <- DNAStringSet()
  F5.Reads.Aligned.Duplex <- DNAStringSet()
  F6.Reads.Aligned.Duplex <- DNAStringSet()
  
  F1.Reads.Consensus.Single <-DNAStringSet()
  F2.Reads.Consensus.Single <-DNAStringSet()
  F3.Reads.Consensus.Single <-DNAStringSet()
  F4.Reads.Consensus.Single <-DNAStringSet()
  F5.Reads.Consensus.Single <-DNAStringSet()
  F6.Reads.Consensus.Single <-DNAStringSet()
  
  F1.Reads.Consensus.Duplex <-DNAStringSet()
  F2.Reads.Consensus.Duplex <-DNAStringSet()
  F3.Reads.Consensus.Duplex <-DNAStringSet()
  F4.Reads.Consensus.Duplex <-DNAStringSet()
  F5.Reads.Consensus.Duplex <-DNAStringSet()
  F6.Reads.Consensus.Duplex <-DNAStringSet()
  

  #Split reads up into groups based on 5' primer (as long as input FASTA file is not empty)
  #Reads are removed from dataset once a primer is found, and only the leftover reads are passed to the next round of primer searching
  #Primers are searched for in the order F1, F2, F3, F4, F5, F6, R1, R2, R3, R4, R5, R6
  if(length(Input) != 0){
    for (f in 1:length(Input)){
      if (vcountPattern(PrimerF1,subseq(Input[f],start=1,end=Search.Distance),fixed=F) == 1){
        F1.Holding <- c(Input[f],F1.Holding)
      }else{
        F1.Leftovers <- c(Input[f],F1.Leftovers)
      }
    }
  }
  if(length(F1.Leftovers) !=0){
    for (f in 1:length(F1.Leftovers)){
      if (vcountPattern(PrimerF2,subseq(F1.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        F2.Holding <- c(F1.Leftovers[f],F2.Holding)
      }else{
        F2.Leftovers <- c(F1.Leftovers[f],F2.Leftovers)
      }
    }
  }
  if(length(F2.Leftovers) !=0){
    for (f in 1:length(F2.Leftovers)){
      if (vcountPattern(PrimerF3,subseq(F2.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        F3.Holding <- c(F2.Leftovers[f],F3.Holding)
      }else{
        F3.Leftovers <- c(F2.Leftovers[f],F3.Leftovers)
      }
    }
  }
  if(length(F3.Leftovers) !=0){
    for (f in 1:length(F3.Leftovers)){
      if (vcountPattern(PrimerF4,subseq(F3.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        F4.Holding <- c(F3.Leftovers[f],F4.Holding)
      }else{
        F4.Leftovers <- c(F3.Leftovers[f],F4.Leftovers)
      }
    }
  }
  if(length(F4.Leftovers) !=0){
    for (f in 1:length(F4.Leftovers)){
      if (vcountPattern(PrimerF5,subseq(F4.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        F5.Holding <- c(F4.Leftovers[f],F5.Holding)
      }else{
        F5.Leftovers <- c(F4.Leftovers[f],F5.Leftovers)
      }
    }
  }
  if(length(F5.Leftovers) !=0){
    for (f in 1:length(F5.Leftovers)){
      if (vcountPattern(PrimerF6,subseq(F5.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        F6.Holding <- c(F5.Leftovers[f],F6.Holding)
      }else{
        F6.Leftovers <- c(F5.Leftovers[f],F6.Leftovers)
      }
    }
  }
  if(length(F6.Leftovers) !=0){
    for (f in 1:length(F6.Leftovers)){
      if (vcountPattern(PrimerR1,subseq(F6.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R1.Holding <- c(F6.Leftovers[f],R1.Holding)
      }else{
        R1.Leftovers <- c(F6.Leftovers[f],R1.Leftovers)
      }
    }
  }
  if(length(R1.Leftovers) !=0){
    for (f in 1:length(R1.Leftovers)){
      if (vcountPattern(PrimerR2,subseq(R1.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R2.Holding <- c(R1.Leftovers[f],R2.Holding)
      }else{
        R2.Leftovers <- c(R1.Leftovers[f],R2.Leftovers)
      }
    }
  }
  if(length(R2.Leftovers) !=0){
    for (f in 1:length(R2.Leftovers)){
      if (vcountPattern(PrimerR3,subseq(R2.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R3.Holding <- c(R2.Leftovers[f],R3.Holding)
      }else{
        R3.Leftovers <- c(R2.Leftovers[f],R3.Leftovers)
      }
    }
  }
  if(length(R3.Leftovers) !=0){
    for (f in 1:length(R3.Leftovers)){
      if (vcountPattern(PrimerR4,subseq(R3.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R4.Holding <- c(R3.Leftovers[f],R4.Holding)
      }else{
        R4.Leftovers <- c(R3.Leftovers[f],R4.Leftovers)
      }
    }
  }
  if(length(R4.Leftovers) !=0){
    for (f in 1:length(R4.Leftovers)){
      if (vcountPattern(PrimerR5,subseq(R4.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R5.Holding <- c(R4.Leftovers[f],R5.Holding)
      }else{
        R5.Leftovers <- c(R4.Leftovers[f],R5.Leftovers)
      }
    }
  }
  if(length(R5.Leftovers) !=0){
    for (f in 1:length(R5.Leftovers)){
      if (vcountPattern(PrimerR6,subseq(R5.Leftovers[f],start=1,end=Search.Distance),fixed=F) == 1){
        R6.Holding <- c(R5.Leftovers[f],R6.Holding)
      }else{
        R6.Leftovers <- c(R5.Leftovers[f],R6.Leftovers)
      }
    }
  }
  
  #reverse compliment reads that contained reverse primers and combine with forward reads
  F1.Reads <- c(F1.Holding,reverseComplement(R1.Holding))
  F2.Reads <- c(F2.Holding,reverseComplement(R2.Holding))
  F3.Reads <- c(F3.Holding,reverseComplement(R3.Holding))
  F4.Reads <- c(F4.Holding,reverseComplement(R4.Holding))
  F5.Reads <- c(F5.Holding,reverseComplement(R5.Holding))
  F6.Reads <- c(F6.Holding,reverseComplement(R6.Holding))
    
  #Separate out singleton reads (e.g. F1-R1) and filter out reads that are too short or too long (expected amplicon length +/- 20 bp)
  F1.Reads.Single <- F1.Reads[which(width(F1.Reads)>249 & width(F1.Reads)<289)] # expected = 269 bp
  F2.Reads.Single <- F2.Reads[which(width(F2.Reads)>224 & width(F2.Reads)<264)] # expected = 244 bp
  F3.Reads.Single <- F3.Reads[which(width(F3.Reads)>234 & width(F3.Reads)<274)] # expected = 254 bp
  F4.Reads.Single <- F4.Reads[which(width(F4.Reads)>258 & width(F4.Reads)<298)] # expected = 278 bp
  F5.Reads.Single <- F5.Reads[which(width(F5.Reads)>233 & width(F5.Reads)<273)] # expected = 253 bp
  F6.Reads.Single <- F6.Reads[which(width(F6.Reads)>240 & width(F6.Reads)<280)] # expected = 260 bp
  
  #Separate out duplex reads (e.g. F1-R2) (expected amplicon +/- 10 bp)
  F1.Reads.Duplex <- F1.Reads[which(width(F1.Reads)>343 & width(F1.Reads)<383)] # expected = 363 bp
  F2.Reads.Duplex <- F2.Reads[which(width(F2.Reads)>311 & width(F2.Reads)<351)] # expected = 331 bp
  F3.Reads.Duplex <- F3.Reads[which(width(F3.Reads)>339 & width(F3.Reads)<379)] # expected = 359 bp
  F4.Reads.Duplex <- F4.Reads[which(width(F4.Reads)>378 & width(F4.Reads)<418)] # expected = 398 bp
  F5.Reads.Duplex <- F5.Reads[which(width(F5.Reads)>341 & width(F5.Reads)<381)] # expected = 361 bp

  #clip adapters and primers from singleton reads
  if (length(F1.Reads.Single)!=0){
    for (j in 1:length(F1.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF1,subseq(F1.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR1.rc,subseq(F1.Reads.Single[j],start=width(F1.Reads.Single[j])-Search.Distance,end=width(F1.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos))) - 1 + width(F1.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F1.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F1.Reads.Single[j],start = start,end = width(F1.Reads.Single[j]) - 46 - width(PrimerR1) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F1.Reads.Single[j],start = 30 + width(PrimerF1) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F1.Reads.Single[j],start = 30 + width(PrimerF1) + 1,end = width(F1.Reads.Single[j]) - 46 - width(PrimerR1) - 1)
        } 
      }
      F1.Reads.Trimmed.Single <- c(F1.Reads.Trimmed.Single,tmp.trimmed)
    }
  }

  if (length(F2.Reads.Single)!=0){
    for (j in 1:length(F2.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF2,subseq(F2.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR2.rc,subseq(F2.Reads.Single[j],start=width(F2.Reads.Single[j])-Search.Distance,end=width(F2.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F2.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F2.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F2.Reads.Single[j],start = start,end = width(F2.Reads.Single[j]) - 46 - width(PrimerR2) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F2.Reads.Single[j],start = 30 + width(PrimerF2) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F2.Reads.Single[j],start = 30 + width(PrimerF2) + 1,end = width(F2.Reads.Single[j]) - 46 - width(PrimerR2) - 1)
        } 
      }
      F2.Reads.Trimmed.Single <- c(F2.Reads.Trimmed.Single,tmp.trimmed)
    }
  }
  
  if (length(F3.Reads.Single)!=0){
    for (j in 1:length(F3.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF3,subseq(F3.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR3.rc,subseq(F3.Reads.Single[j],start=width(F3.Reads.Single[j])-Search.Distance,end=width(F3.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F3.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F3.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F3.Reads.Single[j],start = start,end = width(F3.Reads.Single[j]) - 46 - width(PrimerR3) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F3.Reads.Single[j],start = 30 + width(PrimerF3) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F3.Reads.Single[j],start = 30 + width(PrimerF3) + 1,end = width(F3.Reads.Single[j]) - 46 - width(PrimerR3) - 1)
        }
      }
      F3.Reads.Trimmed.Single <- c(F3.Reads.Trimmed.Single,tmp.trimmed)
    }
  }
  
  if (length(F4.Reads.Single)!=0){
    for (j in 1:length(F4.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF4,subseq(F4.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1 
      rev_pos <- vmatchPattern(PrimerR4.rc,subseq(F4.Reads.Single[j],start=width(F4.Reads.Single[j])-Search.Distance,end=width(F4.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F4.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F4.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F4.Reads.Single[j],start = start,end = width(F4.Reads.Single[j]) - 46 - width(PrimerR4) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F4.Reads.Single[j],start = 30 + width(PrimerF4) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F4.Reads.Single[j],start = 30 + width(PrimerF4) + 1,end = width(F4.Reads.Single[j]) - 46 - width(PrimerR4) - 1)
        }
      }
      F4.Reads.Trimmed.Single <- c(F4.Reads.Trimmed.Single,tmp.trimmed)
    }
  }
  
  if (length(F5.Reads.Single)!=0){
    for (j in 1:length(F5.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF5,subseq(F5.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR5.rc,subseq(F5.Reads.Single[j],start=width(F5.Reads.Single[j])-Search.Distance,end=width(F5.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F5.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F5.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F5.Reads.Single[j],start = start,end = width(F5.Reads.Single[j]) - 46 - width(PrimerR5) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F5.Reads.Single[j],start = 30 + width(PrimerF5) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F5.Reads.Single[j],start = 30 + width(PrimerF5) + 1,end = width(F5.Reads.Single[j]) - 46 - width(PrimerR5) - 1)
        }
      }
      F5.Reads.Trimmed.Single <- c(F5.Reads.Trimmed.Single,tmp.trimmed)
    }
  }
  
  if (length(F6.Reads.Single)!=0){
    for (j in 1:length(F6.Reads.Single)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF6,subseq(F6.Reads.Single[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1 
      rev_pos <- vmatchPattern(PrimerR6.rc,subseq(F6.Reads.Single[j],start=width(F6.Reads.Single[j])-Search.Distance,end=width(F6.Reads.Single[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F6.Reads.Single[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F6.Reads.Single[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F6.Reads.Single[j],start = start,end = width(F6.Reads.Single[j]) - 46 - width(PrimerR6) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F6.Reads.Single[j],start = 30 + width(PrimerF6) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F6.Reads.Single[j],start = 30 + width(PrimerF6) + 1,end = width(F6.Reads.Single[j]) - 46 - width(PrimerR6) - 1)
        } 
      }
      F6.Reads.Trimmed.Single <- c(F6.Reads.Trimmed.Single,tmp.trimmed)
    }
  }

  #clip adapters and primers from duplex reads
  if (length(F1.Reads.Duplex)!=0){
    for (j in 1:length(F1.Reads.Duplex)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF1,subseq(F1.Reads.Duplex[j],start=1,end=Search.Distance),fixed=F) 
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR2.rc,subseq(F1.Reads.Duplex[j],start=width(F1.Reads.Duplex[j])-Search.Distance,end=width(F1.Reads.Duplex[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F1.Reads.Duplex[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F1.Reads.Duplex[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F1.Reads.Duplex[j],start = start,end = width(F1.Reads.Duplex[j]) - 46 - width(PrimerR2) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F1.Reads.Duplex[j],start =+ 30 + width(PrimerF1) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F1.Reads.Duplex[j],start = 30 + width(PrimerF1) + 1,end = width(F1.Reads.Duplex[j]) - 46 - width(PrimerR2) - 1)
        } 
      }
      F1.Reads.Trimmed.Duplex <- c(F1.Reads.Trimmed.Duplex,tmp.trimmed)
    }
  }
  
  if (length(F2.Reads.Duplex)!=0){
    for (j in 1:length(F2.Reads.Duplex)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF2,subseq(F2.Reads.Duplex[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR3.rc,subseq(F2.Reads.Duplex[j],start=width(F2.Reads.Duplex[j])-Search.Distance,end=width(F2.Reads.Duplex[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F2.Reads.Duplex[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F2.Reads.Duplex[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F2.Reads.Duplex[j],start = start,end = width(F2.Reads.Duplex[j]) - 46 - width(PrimerR3) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F2.Reads.Duplex[j],start = 30 + width(PrimerF2) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F2.Reads.Duplex[j],start = 30 + width(PrimerF2) + 1,end = width(F2.Reads.Duplex[j]) - 46 - width(PrimerR3) - 1)
        } 
      }
      F2.Reads.Trimmed.Duplex <- c(F2.Reads.Trimmed.Duplex,tmp.trimmed)
    }
  }
  
  if (length(F3.Reads.Duplex)!=0){
    for (j in 1:length(F3.Reads.Duplex)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF3,subseq(F3.Reads.Duplex[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR4.rc,subseq(F3.Reads.Duplex[j],start=width(F3.Reads.Duplex[j])-Search.Distance,end=width(F3.Reads.Duplex[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F3.Reads.Duplex[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F3.Reads.Duplex[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F3.Reads.Duplex[j],start = start,end = width(F3.Reads.Duplex[j]) - 46 - width(PrimerR4) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F3.Reads.Duplex[j],start = 30 + width(PrimerF3) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F3.Reads.Duplex[j],start = 30 + width(PrimerF3) + 1,end = width(F3.Reads.Duplex[j]) - 46 - width(PrimerR4) - 1)
        } 
      }
      F3.Reads.Trimmed.Duplex <- c(F3.Reads.Trimmed.Duplex,tmp.trimmed)
    }
  }
  
  if (length(F4.Reads.Duplex)!=0){
    for (j in 1:length(F4.Reads.Duplex)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF4,subseq(F4.Reads.Duplex[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1
      rev_pos <- vmatchPattern(PrimerR5.rc,subseq(F4.Reads.Duplex[j],start=width(F4.Reads.Duplex[j])-Search.Distance,end=width(F4.Reads.Duplex[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F4.Reads.Duplex[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F4.Reads.Duplex[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F4.Reads.Duplex[j],start = start,end = width(F4.Reads.Duplex[j]) - 46 - width(PrimerR5) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F4.Reads.Duplex[j],start = 30 + width(PrimerF4) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F4.Reads.Duplex[j],start = 30 + width(PrimerF4) + 1,end = width(F4.Reads.Duplex[j]) - 46 - width(PrimerR5) - 1)
        } 
      }
      F4.Reads.Trimmed.Duplex <- c(F4.Reads.Trimmed.Duplex,tmp.trimmed)
    }
  }
  
  if (length(F5.Reads.Duplex)!=0){
    for (j in 1:length(F5.Reads.Duplex)){
      tmp.trimmed <- DNAStringSet()
      fwd_pos <- vmatchPattern(PrimerF5,subseq(F5.Reads.Duplex[j],start=1,end=Search.Distance),fixed=F)
      start <- as.numeric(max(end(fwd_pos)))+1 
      rev_pos <- vmatchPattern(PrimerR6.rc,subseq(F5.Reads.Duplex[j],start=width(F5.Reads.Duplex[j])-Search.Distance,end=width(F5.Reads.Duplex[j])),fixed=F)
      end <- as.numeric(min(start(rev_pos)))-1 + width(F5.Reads.Duplex[j]) - Search.Distance
      if (end - start > 1){
        if (is.infinite(start)==FALSE & is.infinite(end)==FALSE){ #both primers found
          tmp.trimmed <- subseq(F5.Reads.Duplex[j],start = start,end = end -1)
        } else if (is.infinite(start)==FALSE & is.infinite(end)==TRUE){ #forward found, reverse not found
          tmp.trimmed <- subseq(F5.Reads.Duplex[j],start = start,end = width(F5.Reads.Duplex[j]) - 46 - width(PrimerR6) - 1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==FALSE){ #forward not found, reverse found
          tmp.trimmed <- subseq(F5.Reads.Duplex[j],start = 30 + width(PrimerF5) + 1,end = end -1)
        } else if (is.infinite(start)==TRUE & is.infinite(end)==TRUE){ #neither primer found
          tmp.trimmed <- subseq(F5.Reads.Duplex[j],start = 30 + width(PrimerF5) + 1,end = width(F5.Reads.Duplex[j]) - 46 - width(PrimerR6) - 1)
        } 
      }
      F5.Reads.Trimmed.Duplex <- c(F5.Reads.Trimmed.Duplex,tmp.trimmed)
    }
  }
  
  #Filter out reads greater or less than target amplicon by 2bp
  F1.Reads.Trimmed.Single <- F1.Reads.Trimmed.Single[which(width(F1.Reads.Trimmed.Single) >= 142 & width(F1.Reads.Trimmed.Single) <= 146)] # expected = 144 bp
  F2.Reads.Trimmed.Single <- F2.Reads.Trimmed.Single[which(width(F2.Reads.Trimmed.Single) >= 118 & width(F2.Reads.Trimmed.Single) <= 122)] # expected = 120 bp
  F3.Reads.Trimmed.Single <- F3.Reads.Trimmed.Single[which(width(F3.Reads.Trimmed.Single) >= 128 & width(F3.Reads.Trimmed.Single) <= 132)] # expected = 130 bp
  F4.Reads.Trimmed.Single <- F4.Reads.Trimmed.Single[which(width(F4.Reads.Trimmed.Single) >= 155 & width(F4.Reads.Trimmed.Single) <= 159)] # expected = 157 bp
  F5.Reads.Trimmed.Single <- F5.Reads.Trimmed.Single[which(width(F5.Reads.Trimmed.Single) >= 134 & width(F5.Reads.Trimmed.Single) <= 138)] # expected = 136 bp
  F6.Reads.Trimmed.Single <- F6.Reads.Trimmed.Single[which(width(F6.Reads.Trimmed.Single) >= 131 & width(F6.Reads.Trimmed.Single) <= 135)] # expected = 133 bp
  
  F1.Reads.Trimmed.Duplex <- F1.Reads.Trimmed.Duplex[which(width(F1.Reads.Trimmed.Duplex) >= 236 & width(F1.Reads.Trimmed.Duplex) <= 240)] # expected = 238 bp
  F2.Reads.Trimmed.Duplex <- F2.Reads.Trimmed.Duplex[which(width(F2.Reads.Trimmed.Duplex) >= 205 & width(F2.Reads.Trimmed.Duplex) <= 209)] # expected = 207 bp
  F3.Reads.Trimmed.Duplex <- F3.Reads.Trimmed.Duplex[which(width(F3.Reads.Trimmed.Duplex) >= 233 & width(F3.Reads.Trimmed.Duplex) <= 237)] # expected = 235 bp
  F4.Reads.Trimmed.Duplex <- F4.Reads.Trimmed.Duplex[which(width(F4.Reads.Trimmed.Duplex) >= 275 & width(F4.Reads.Trimmed.Duplex) <= 279)] # expected = 277 bp
  F5.Reads.Trimmed.Duplex <- F5.Reads.Trimmed.Duplex[which(width(F5.Reads.Trimmed.Duplex) >= 242 & width(F5.Reads.Trimmed.Duplex) <= 246)] # expected = 244 bp

  #Align reads within each fragment and take consensus, then de-gap consensus and replicate it to match number of origianl reads (provides "weighted" consensus seqs)
  #F1 singletons
  if(length(F1.Reads.Trimmed.Single) >=3){
    F1.Reads.Aligned.Single <- DNAStringSet(muscle(F1.Reads.Trimmed.Single))
    F1.Reads.Consensus.Single <- DNAStringSet(consensusString(F1.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F1.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F1.Reads.Consensus.Single))
    F1.Reads.Consensus.Single <- DNAStringSet(rep(F1.Reads.Consensus.Single,times=length(F1.Reads.Aligned.Single)))
  }

  #F1 duplex
  if(length(F1.Reads.Trimmed.Duplex) >=3){
    F1.Reads.Aligned.Duplex <- DNAStringSet(muscle(F1.Reads.Trimmed.Duplex))
    F1.Reads.Consensus.Duplex <- DNAStringSet(consensusString(F1.Reads.Aligned.Duplex, ambiguityMap = "N", threshold = 0.5))
    F1.Reads.Consensus.Duplex <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F1.Reads.Consensus.Duplex))
    F1.Reads.Consensus.Duplex <- DNAStringSet(rep(F1.Reads.Consensus.Duplex,times=length(F1.Reads.Aligned.Duplex)))
  }

  #F2 singletons
  if(length(F2.Reads.Trimmed.Single) >=3){
    F2.Reads.Aligned.Single <- DNAStringSet(muscle(F2.Reads.Trimmed.Single))
    F2.Reads.Consensus.Single <- DNAStringSet(consensusString(F2.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F2.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F2.Reads.Consensus.Single))
    F2.Reads.Consensus.Single <- DNAStringSet(rep(F2.Reads.Consensus.Single,times=length(F2.Reads.Aligned.Single)))
  }

  #F2 duplex
  if(length(F2.Reads.Trimmed.Duplex) >=3){
    F2.Reads.Aligned.Duplex <- DNAStringSet(muscle(F2.Reads.Trimmed.Duplex))
    F2.Reads.Consensus.Duplex <- DNAStringSet(consensusString(F2.Reads.Aligned.Duplex, ambiguityMap = "N", threshold = 0.5))
    F2.Reads.Consensus.Duplex <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F2.Reads.Consensus.Duplex))
    F2.Reads.Consensus.Duplex <- DNAStringSet(rep(F2.Reads.Consensus.Duplex,times=length(F2.Reads.Aligned.Duplex)))
  }

  #F3 singletons
  if(length(F3.Reads.Trimmed.Single) >=3){
    F3.Reads.Aligned.Single <- DNAStringSet(muscle(F3.Reads.Trimmed.Single))
    F3.Reads.Consensus.Single <- DNAStringSet(consensusString(F3.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F3.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F3.Reads.Consensus.Single))
    F3.Reads.Consensus.Single <- DNAStringSet(rep(F3.Reads.Consensus.Single,times=length(F3.Reads.Aligned.Single)))
  }

  #F3 duplex
  if(length(F3.Reads.Trimmed.Duplex) >=3){
    F3.Reads.Aligned.Duplex <- DNAStringSet(muscle(F3.Reads.Trimmed.Duplex))
    F3.Reads.Consensus.Duplex <- DNAStringSet(consensusString(F3.Reads.Aligned.Duplex, ambiguityMap = "N", threshold = 0.5))
    F3.Reads.Consensus.Duplex <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F3.Reads.Consensus.Duplex))
    F3.Reads.Consensus.Duplex <- DNAStringSet(rep(F3.Reads.Consensus.Duplex,times=length(F3.Reads.Aligned.Duplex)))
  }

  #F4 singletons
  if(length(F4.Reads.Trimmed.Single) >=3){
    F4.Reads.Aligned.Single <- DNAStringSet(muscle(F4.Reads.Trimmed.Single))
    F4.Reads.Consensus.Single <- DNAStringSet(consensusString(F4.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F4.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F4.Reads.Consensus.Single))
    F4.Reads.Consensus.Single <- DNAStringSet(rep(F4.Reads.Consensus.Single,times=length(F4.Reads.Aligned.Single)))
  }

  #F4 duplex
  if(length(F4.Reads.Trimmed.Duplex) >=3){
    F4.Reads.Aligned.Duplex <- DNAStringSet(muscle(F4.Reads.Trimmed.Duplex))
    F4.Reads.Consensus.Duplex <- DNAStringSet(consensusString(F4.Reads.Aligned.Duplex, ambiguityMap = "N", threshold = 0.5))
    F4.Reads.Consensus.Duplex <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F4.Reads.Consensus.Duplex))
    F4.Reads.Consensus.Duplex <- DNAStringSet(rep(F4.Reads.Consensus.Duplex,times=length(F4.Reads.Aligned.Duplex)))
  }

  #F5 singletons
  if(length(F5.Reads.Trimmed.Single) >=3){
    F5.Reads.Aligned.Single <- DNAStringSet(muscle(F5.Reads.Trimmed.Single))
    F5.Reads.Consensus.Single <- DNAStringSet(consensusString(F5.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F5.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F5.Reads.Consensus.Single))
    F5.Reads.Consensus.Single <- DNAStringSet(rep(F5.Reads.Consensus.Single,times=length(F5.Reads.Aligned.Single)))
  }

  #F5 duplex
  if(length(F5.Reads.Trimmed.Duplex) >=3){
    F5.Reads.Aligned.Duplex <- DNAStringSet(muscle(F5.Reads.Trimmed.Duplex))
    F5.Reads.Consensus.Duplex <- DNAStringSet(consensusString(F5.Reads.Aligned.Duplex, ambiguityMap = "N", threshold = 0.5))
    F5.Reads.Consensus.Duplex <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F5.Reads.Consensus.Duplex))
    F5.Reads.Consensus.Duplex <- DNAStringSet(rep(F5.Reads.Consensus.Duplex,times=length(F5.Reads.Aligned.Duplex)))
  }

  #F6 singletons
  if(length(F6.Reads.Trimmed.Single) >=3){
    F6.Reads.Aligned.Single <- DNAStringSet(muscle(F6.Reads.Trimmed.Single))
    F6.Reads.Consensus.Single <- DNAStringSet(consensusString(F6.Reads.Aligned.Single, ambiguityMap = "N", threshold = 0.5))
    F6.Reads.Consensus.Single <- DNAStringSet(gsub(pattern = "-", replacement = "", x = F6.Reads.Consensus.Single))
    F6.Reads.Consensus.Single <- DNAStringSet(rep(F6.Reads.Consensus.Single,times=length(F6.Reads.Aligned.Single)))
  }

  #Remove consensus sequences with more than 2 N's
  if(length(F1.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F1.Reads.Consensus.Single[1], "N")) > 2){
      F1.Reads.Consensus.Single <- DNAStringSet()
    }
  }
  
  if(length(F2.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F2.Reads.Consensus.Single[1], "N")) > 2){
      F2.Reads.Consensus.Single <- DNAStringSet()
    }
  }

  if(length(F3.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F3.Reads.Consensus.Single[1], "N")) > 2){
      F3.Reads.Consensus.Single <- DNAStringSet()
    }
  }
  
  if(length(F4.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F4.Reads.Consensus.Single[1], "N")) > 2){
      F4.Reads.Consensus.Single <- DNAStringSet()
    }
  }
  
  if(length(F5.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F5.Reads.Consensus.Single[1], "N")) > 2){
      F5.Reads.Consensus.Single <- DNAStringSet()
    }
  }
  
  if(length(F6.Reads.Consensus.Single)!=0){
    if(as.numeric(letterFrequency(F6.Reads.Consensus.Single[1], "N")) > 2){
      F6.Reads.Consensus.Single <- DNAStringSet()
    }
  }
  
  if(length(F1.Reads.Consensus.Duplex)!=0){
    if(as.numeric(letterFrequency(F1.Reads.Consensus.Duplex[1], "N")) > 2){
      F1.Reads.Consensus.Duplex <- DNAStringSet()
    }
  }
  
  if(length(F2.Reads.Consensus.Duplex)!=0){
    if(as.numeric(letterFrequency(F2.Reads.Consensus.Duplex[1], "N")) > 2){
      F2.Reads.Consensus.Duplex <- DNAStringSet()
    }
  }
  
  if(length(F3.Reads.Consensus.Duplex)!=0){
    if(as.numeric(letterFrequency(F3.Reads.Consensus.Duplex[1], "N")) > 2){
      F3.Reads.Consensus.Duplex <- DNAStringSet()
    }
  }
  
  if(length(F4.Reads.Consensus.Duplex)!=0){
    if(as.numeric(letterFrequency(F4.Reads.Consensus.Duplex[1], "N")) > 2){
      F4.Reads.Consensus.Duplex <- DNAStringSet()
    }
  }
  
  if(length(F5.Reads.Consensus.Duplex)!=0){
    if(as.numeric(letterFrequency(F5.Reads.Consensus.Duplex[1], "N")) > 2){
      F5.Reads.Consensus.Duplex <- DNAStringSet()
    }
  }

  #Add "+" to the beginning of each fragment's consensus sequence to force them into alignment
  F1.Final.Single <- F1.Reads.Consensus.Single
  if(length(F1.Final.Single)!=0){names(F1.Final.Single) <- paste("F1_Single_",seq(from = 1, to = length(F1.Final.Single), by = 1),sep = "")}
  
  F2.Leader.Single <- DNAStringSet(rep(paste(rep("N",118),collapse=""),length(F2.Reads.Consensus.Single)))
  F2.Final.Single <- DNAStringSet(Map(c,F2.Leader.Single,F2.Reads.Consensus.Single))
  if(length(F2.Final.Single)!=0){names(F2.Final.Single) <- paste("F2_Single_",seq(from = 1, to = length(F2.Final.Single), by = 1),sep = "")}
  
  F3.Leader.Single <- DNAStringSet(rep(paste(rep("N",195),collapse=""),length(F3.Reads.Consensus.Single)))
  F3.Final.Single <- DNAStringSet(Map(c,F3.Leader.Single,F3.Reads.Consensus.Single))
  if(length(F3.Final.Single)!=0){names(F3.Final.Single) <- paste("F3_Single_",seq(from = 1, to = length(F3.Final.Single), by = 1),sep = "")}
  
  F4.Leader.Single <- DNAStringSet(rep(paste(rep("N",273),collapse=""),length(F4.Reads.Consensus.Single)))
  F4.Final.Single <- DNAStringSet(Map(c,F4.Leader.Single,F4.Reads.Consensus.Single))
  if(length(F4.Final.Single)!=0){names(F4.Final.Single) <- paste("F4_Single_",seq(from = 1, to = length(F4.Final.Single), by = 1),sep = "")}
  
  F5.Leader.Single <- DNAStringSet(rep(paste(rep("N",414),collapse=""),length(F5.Reads.Consensus.Single)))
  F5.Final.Single <- DNAStringSet(Map(c,F5.Leader.Single,F5.Reads.Consensus.Single))
  if(length(F5.Final.Single)!=0){names(F5.Final.Single) <-paste("F5_Single_",seq(from = 1, to = length(F5.Final.Single), by = 1),sep = "")}
  
  F6.Leader.Single <- DNAStringSet(rep(paste(rep("N",525),collapse=""),length(F6.Reads.Consensus.Single)))
  F6.Final.Single <- DNAStringSet(Map(c,F6.Leader.Single,F6.Reads.Consensus.Single))
  if(length(F6.Final.Single)!=0){names(F6.Final.Single) <- paste("F6_Single_",seq(from = 1, to = length(F6.Final.Single), by = 1),sep = "")}
  
  F1.Final.Duplex <- F1.Reads.Consensus.Duplex
  if(length(F1.Final.Duplex)!=0){names(F1.Final.Duplex) <- paste("F1_Duplex_",seq(from = 1, to = length(F1.Final.Duplex), by = 1),sep = "")}
  
  F2.Leader.Duplex <- DNAStringSet(rep(paste(rep("N",118),collapse=""),length(F2.Reads.Consensus.Duplex)))
  F2.Final.Duplex <- DNAStringSet(Map(c,F2.Leader.Duplex,F2.Reads.Consensus.Duplex))
  if(length(F2.Final.Duplex)!=0){names(F2.Final.Duplex) <- paste("F2_Duplex_",seq(from = 1, to = length(F2.Final.Duplex), by = 1),sep = "")}
  
  F3.Leader.Duplex <- DNAStringSet(rep(paste(rep("N",195),collapse=""),length(F3.Reads.Consensus.Duplex)))
  F3.Final.Duplex <- DNAStringSet(Map(c,F3.Leader.Duplex,F3.Reads.Consensus.Duplex))
  if(length(F3.Final.Duplex)!=0){names(F3.Final.Duplex) <- paste("F3_Duplex_",seq(from = 1, to = length(F3.Final.Duplex), by = 1),sep = "")}
  
  F4.Leader.Duplex <- DNAStringSet(rep(paste(rep("N",273),collapse=""),length(F4.Reads.Consensus.Duplex)))
  F4.Final.Duplex <- DNAStringSet(Map(c,F4.Leader.Duplex,F4.Reads.Consensus.Duplex))
  if(length(F4.Final.Duplex)!=0){names(F4.Final.Duplex) <- paste("F4_Duplex_",seq(from = 1, to = length(F4.Final.Duplex), by = 1),sep = "")}
  
  F5.Leader.Duplex <- DNAStringSet(rep(paste(rep("N",414),collapse=""),length(F5.Reads.Consensus.Duplex)))
  F5.Final.Duplex <- DNAStringSet(Map(c,F5.Leader.Duplex,F5.Reads.Consensus.Duplex))
  if(length(F5.Final.Duplex)!=0){names(F5.Final.Duplex) <-paste("F5_Duplex_",seq(from = 1, to = length(F5.Final.Duplex), by = 1),sep = "")}
  
  #combine all reads into a single file
  All.Final <- DNAStringSet(c(F1.Final.Single,F1.Final.Duplex,F2.Final.Single,F2.Final.Duplex,F3.Final.Single,F3.Final.Duplex,F4.Final.Single,F4.Final.Duplex,F5.Final.Single,F5.Final.Duplex,F6.Final.Single))

  #generate consensus sequence of reads
  if(length(All.Final) !=0){
    All.Final <- DNAStringSet(gsub(pattern = "N", replacement = "+", x = All.Final))
    Consensus <- DNAStringSet(ConsensusSequence(All.Final, ignoreNonBases = TRUE, noConsensusChar = "N",threshold = 0.6, minInformation = 0.4))
    Consensus <- DNAStringSet(gsub(pattern = "-", replacement = "N", x = Consensus))
    Sample.Taxonomy <- as.character(Sample.Taxonomy.Table[grep(pattern = Sample.Name, x = Sample.Taxonomy.Table$ProcessID),4])
    names(Consensus) <- paste(Sample.Name,Sample.Taxonomy,sep = "|")
    Reads.Final <- DNAStringSet(c(Consensus,All.Final))
    writeXStringSet(Reads.Final,sprintf("%s_ASSEMBLED.fasta",Sample.Name))

    #add consensus sequence to final output file
    output <- c(Consensus,output)
  }
  
  #Track individual reads per fragment
  F1.Reads.Fragment <- length(c(F1.Final.Single,F1.Final.Duplex))
  F2.Reads.Fragment <- length(c(F2.Final.Single,F2.Final.Duplex,F1.Final.Duplex))
  F3.Reads.Fragment <- length(c(F3.Final.Single,F3.Final.Duplex,F2.Final.Duplex))
  F4.Reads.Fragment <- length(c(F4.Final.Single,F4.Final.Duplex,F3.Final.Duplex))
  F5.Reads.Fragment <- length(c(F5.Final.Single,F5.Final.Duplex,F4.Final.Duplex))
  F6.Reads.Fragment <- length(c(F6.Final.Single,F5.Final.Duplex))
  tmp <- data.frame("Sample" = Sample.Name,"F1(1-144)" = F1.Reads.Fragment,"F2(119-238)" = F2.Reads.Fragment,"F3(196-325)" = F3.Reads.Fragment,"F4(274-430)" = F4.Reads.Fragment,"F5(415-550)" = F5.Reads.Fragment,"F6(526-658)" = F6.Reads.Fragment)
  tmp.duplex <- data.frame("Sample" = Sample.Name,"F1(1-144)" = length(F1.Final.Duplex),"F2(119-238)" = length(F2.Final.Duplex),"F3(196-325)" = length(F3.Final.Duplex),"F4(274-430)" = length(F4.Final.Duplex),"F5(415-550)" = length(F5.Final.Duplex),"F6: 526-658" = 0)
  Reads.Table <- rbind(Reads.Table,tmp)
  Duplex.Table <- rbind(Duplex.Table,tmp.duplex)
} 

#Write all final outputs
writeXStringSet(output,"Original Output.fas")
write.csv(Reads.Table,"Coverage Per Fragment.csv",row.names=FALSE)
write.csv(Duplex.Table,"Duplex Per Fragment.csv",row.names=FALSE)

  
  
  
