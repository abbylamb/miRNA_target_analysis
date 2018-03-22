#########Playing around with adding new info from flybase (expression,etc) but will need to add FBids

ids_miranda<-subset(all, select = c(gene_id, gene_symbol, transcript_id))
annotations <- read.delim(file = "flybase_annotation_table.txt", sep = "\t")
annotations <- data.frame(lapply(annotations, as.character), stringsAsFactors = FALSE)

##need annotation ID to match flybase annotation table, so I'm taking the "-R_" off each transcript_id
ids_miranda <- mutate(ids_miranda, annotation_id = substr(transcript_id,1,nchar(transcript_id)-3))
head(ids_miranda)
dim(ids_miranda)
class(ids_miranda$transcript_id)
head(annotations)
colnames(annotations) <- c("gene_symbol", "primary_FBgn", "secondary_FBgns", "annotation_ID", "secondary_annotation_IDs") 
annotations_reduced <- subset(annotations, select = c(primary_FBgn, annotation_ID))
head(annotations_reduced)
##want to add FBgn numbers to each id in my database
testing_ljoin <- left_join(ids_miranda, annotations_reduced, by = "annotation_ID")
annotations_reduced <- data.frame(lapply(annotations_reduced, as.character), stringsAsFactors=FALSE)
ids_miranda <- data.frame(lapply(ids_miranda, as.character), stringsAsFactors=FALSE)
class(ids_miranda$annotation_ID)
colnames(ids_miranda) <- c("gene_ID", "gene_symbol", "transcript_ID", "annotation_ID")
write.csv(testing_ljoin, "added_FBgn.csv")
##quite a few lines didn't have matches in annotation ID. 
## need to figure out how to use secondary annotation IDs

###the next 5 lines are garbage
#ljoin2 <- data.frame(testing_ljoin) 
#names(ljoin2) <- c("gene_ID", "gene_symbol", "transcript_ID", "secondary_annotation_IDs", "primary_FBgn")
#ljoin2 <- left_join(ljoin2, annotations, by = "secondary_annotation_IDs")
#head(ljoin2)
#write.csv(ljoin2, "maybe_better_added_FBgn.csv")

install.packages("tidyr")
library("tidyr")

##need to split secondary annotation IDs into individual columns
septest <- data.frame(annotations, stringsAsFactors = FALSE)
split_IDs <- septest %>% separate(secondary_annotation_IDs, c("secondary_ID", "tertiary_ID", "quaternary_ID", "5th_ID", "6th_ID", "7th_ID", "8th_ID", "9th_ID", "10th_ID", "11th_ID", "12th_ID", "13th_ID"), ",")
## Kept adding columns until I stopped getting warnings about additional pieces being discarded.
head(split_IDs)
colnames(split_IDs)
write.csv(split_IDs, "split_IDs_flybase_Feb14.csv")

## So far, "testingljoin" is the best I have as far as a combination of Flybase and miRanda names
## need to find a way to get flybase secondary annotation ids matched up to miRanda list
## in order to match a FBgn number to every gene.
counter <- 1
looptest <- data.frame(testing_ljoin, stringsAsFactors = FALSE)
for (i in looptest$primary_FBgn){
  
  if (is.na(looptest$primary_FBgn[counter]) == TRUE & 
      looptest$annotation_ID[counter]  %in% split_IDs$secondary_ID == TRUE
  )
  {
    
    looptest$primary_FBgn[counter] <- split_IDs$primary_FBgn[split_IDs$secondary_ID == looptest$annotation_ID[counter]]
  }
  
  
  counter <- counter + 1
}
with_secondary <- data.frame(looptest, stringsAsFactors = FALSE)
write.csv(with_secondary, "with_secondary.csv")

### Now that I've added secondary annotation IDs, need to add tertiary annotation IDs
with_tertiary <- data.frame(with_secondary, stringsAsFactors = FALSE)
counter <- 1
for (i in with_tertiary$primary_FBgn){
  
  if (is.na(with_tertiary$primary_FBgn[counter]) == TRUE & with_tertiary$annotation_ID[counter]  %in% split_IDs$tertiary_ID == TRUE)
  {
    with_tertiary$primary_FBgn[counter] <- split_IDs$primary_FBgn[split_IDs$tertiary_ID == with_tertiary$annotation_ID[counter] & is.na(split_IDs$tertiary_ID) == FALSE]
  }
  counter <- counter + 1
}
head(with_tertiary)
write.csv(with_tertiary, "with_tertiary.csv")

## adding quaternary annotation ids. down to only 50 unmatched annotation IDs so far!
with_quaternary <- data.frame(with_tertiary, stringsAsFactors = FALSE)
counter <- 1
for (i in with_quaternary$primary_FBgn){
  
  if (is.na(with_quaternary$primary_FBgn[counter]) == TRUE & with_quaternary$annotation_ID[counter]  %in% split_IDs$quaternary_ID == TRUE)
  {
    with_quaternary$primary_FBgn[counter] <- split_IDs$primary_FBgn[split_IDs$quaternary_ID == with_quaternary$annotation_ID[counter] & is.na(split_IDs$quaternary_ID) == FALSE]
  }
  counter <- counter + 1
}
head(with_quaternary)
write.csv(with_quaternary, "with_quaternary.csv")

## adding 5th annotation ids. down to only 33 unmatched annotation IDs so far!
with_5th <- data.frame(with_quaternary, stringsAsFactors = FALSE)
counter <- 1
for (i in with_5th$primary_FBgn){
  
  if (is.na(with_5th$primary_FBgn[counter]) == TRUE & with_5th$annotation_ID[counter]  %in% split_IDs$`5th_ID` == TRUE)
  {
    with_5th$primary_FBgn[counter] <- split_IDs$primary_FBgn[split_IDs$`5th_ID` == with_5th$annotation_ID[counter] & is.na(split_IDs$`5th_ID`) == FALSE]
  }
  counter <- counter + 1
}
head(with_5th)
write.csv(with_5th, "with_5th.csv")

### 21 unpaired annotation IDs remain. May just fix them manually.
### it turns out most of these 21 aren't even genes anymore. they've been removed, but still have ids

dim(with_5th)
no_repeat_with_5th <- unique(with_5th)
dim(no_repeat_with_5th)
write.csv(no_repeat_with_5th, "FBgns_miranda_no_repeats.csv")
###opened up the csv file above and manually annotated the remaining few fbgns.
###there are still repeats of genes because of different transcript ids, so I'll get rid of those too.

miranda_fbgn <- read.csv("miranda_genes_w_FBgn.csv")
head(miranda_fbgn)
miranda_fbgn <- subset(miranda_fbgn, select = c(gene_ID, primary_FBgn))
head(miranda_fbgn)
dim(miranda_fbgn)
miranda_fbgn <- unique(miranda_fbgn)
dim(miranda_fbgn)
write.csv(miranda_fbgn, "unique_miranda_w_fbgn.csv")
colnames(miranda_fbgn) <- c("gene_id", "fbgn")
### should now only have unique miranda gene IDs and their corresponding FBgn numbers

## add FBgn data to main dataset
all_with_fbgn <- full_join(all, miranda_fbgn, by="gene_id")
fullset_with_fbgn <- full_join(fullset, miranda_fbgn, by="gene_id")
write.csv(fullset_with_fbgn, "fullset_screen_piggenes_fbgn.csv")