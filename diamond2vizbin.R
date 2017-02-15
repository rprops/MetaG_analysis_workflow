### Script for translating Diamond classification output to vizbin annotation
### Output is a file with contig ID and label at the requested taxonomic rank
### The user will still need to sort the list so that the rows exactly match
### with the contigs used to generate the binning.
userprefs <- commandArgs(trailingOnly = TRUE)

library("dplyr")

file <- userprefs[1]
reference <- userprefs[2]

file1 <- read.delim(file,
                       stringsAsFactors = FALSE, header=FALSE)
reference <- read.delim(reference,
                        stringsAsFactors = FALSE, header=FALSE)

tmp <- anti_join(reference, file1, by = "V1")
tmp <- data.frame(V1 = tmp, V2 =rep("UNCLASSIFIED", nrow(reference) - nrow(file1)))
result <- rbind(file1, tmp)
result <- result[match(reference$V1, as.character(result$V1)),]
write.csv(data.frame(label = as.character(result$V2)), file="Annotation_vizbin.csv", row.names = FALSE)