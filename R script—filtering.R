setwd("/mnt/Files_1/sprosser/R_Working_Directory")
library(Biostrings)

#Compare BLAST results to BOLD taxonomy and filter out reads that don't match 

#initialize some variables
Order <- read.table("samples.txt", header = TRUE, sep = "\t")
results <- list()

#consolidate and reformat BLAST results
for(f in list.files(pattern = "*assignments.txt")){
  sample <- gsub("_READS_tax_assignments.txt","",f)
  df <- read.csv(file=f, sep="\t", col.names=c("Read","Identification","E-value","Percent_ID","Coverage","Reference_Process_ID"), header=TRUE)
  if(nrow(df) == 0){
    df <- data.frame("Sample" = sample, "Read" = "failed", "BLAST.Order" = "failed")
    results <- rbind(results,df)
  } else{
    df$Sample <- rep(sample,times=length(df[,1]))
    df$Identification <- gsub("No blast hit", "no match;no match;no match;no match;no match;no match",df$Identification)
    df <- data.frame(df,do.call(rbind, strsplit(as.character(df$Identification),";",fixed=TRUE)))
    if(length(df) < 10){
      df$temp <- ""
    }
    df <- df[,c(7,1,10)]
    colnames(df) <- c("Sample","Read","BLAST.Order")
    results <- rbind(results,df)
    }
}

#compare BLAST results to BOLD ID at ordinal level, output table of results
results$BOLD.Order <- Order[match(results$Sample, Order$ProcessID),3]
results$Red.Flag <- ifelse(as.character(results$BLAST.Order) != as.character(results$BOLD.Order),"Yes","No")
write.table(x = results, file = "BLAST_Filtering_Summary.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#read in FASTA files and filter out reads that don't match at order level
targets <- results[results$Red.Flag=="Yes",c(1,2,5)]
for(i in list.files(pattern = "*.fasta")){
  sample <- gsub("_READS.fasta","",i)
  seqs <- readDNAStringSet(i)
  targets.sample <- targets[targets$Sample==sample,]
  target.indices <- which(!(names(seqs) %in% as.character(targets.sample$Read)))
  nontarget.indices <- which((names(seqs) %in% as.character(targets.sample$Read)))
  output.fasta <- seqs[target.indices]
  output.filtered <- seqs[nontarget.indices]
  writeXStringSet(output.fasta,paste(sample,"_TARGETS.fasta",sep = ""))
  writeXStringSet(output.filtered,paste(sample,"_NONTARGETS.fasta",sep = ""))
}






















