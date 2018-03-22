##Read in both datasets from miRanda and combine "conserved" and "nonconserved" sets

getwd()

cons <- read.delim("miranda_c.txt", sep = "\t")
ncons <- read.delim("miranda_nc.txt", sep = "\t")

library("dplyr")

all <- full_join(cons, ncons)
nrow(all) == nrow(cons)+nrow(ncons)
head(all)

##read in list of all pigment genes, 
##make separate dataframes of genes that darken pigmentation (e.g. yellow) 
##or lighten pigmentation (e.g. ebony)
##updated to include all pigment genes I could find in Flybase that affect adult cuticle pigment
piggenes <- read.csv("gene_list_expanded_dsx_both.csv")
darkgenes <- subset(piggenes, gene_fxn_pigm == 'D' | gene_fxn_pigm == 'B')
lightgenes <- subset(piggenes, gene_fxn_pigm == 'L' | gene_fxn_pigm == 'B')
View(piggenes)


##read in overexpression screen results for miRNAs (the one with names that match miRanda database)
##make separate dataframes of miRs that, when overexpressed, darken or lighten pigmentation
dmemiRs <- read.csv("tested_mirs_match_miranda.csv")
darkmirs <- subset(dmemiRs, OE_pigment_PT == 'B' | OE_pigment_PT == 'D')
lightmirs <- subset(dmemiRs, OE_pigment_PT == 'B' | OE_pigment_PT == 'L')
noeffectmirs<- subset(dmemiRs, OE_pigment_PT == 'N')

## rename columns of miRNA names to match between databases, then check
names(dmemiRs)[1] <- names(all[2])
names(all)
names(dmemiRs)

## combine dataframes of screen results and miRanda data (dmemiRs and all)
all_w_screen <- full_join(all, dmemiRs, by = "mirna_name")
head(all_w_screen)
dim(all_w_screen)
dim(all)
##there's 9 more rows in the joined set. find what's different. 
##but I think setdiff won't work with 2 different sets of columns

##Find all unique miRs from all and from all_w_screen. Compare the 2 lists to see what's different
all_mirs <- unique(all$mirna_name)
my_mirs <- unique(dmemiRs$mirna_name)
all_w_screen_mirs <- unique(all_w_screen$mirna_name)
setdiff(all_w_screen_mirs, all_mirs)
## above is a list of all mirs that are in my screen but not in the miranda database.
## Since correcting the naming issues, now the only mismatches are the miRs that 
## I screened that are not in miRanda (3000s and 4000s). Good!

## For the time being, I'll just go on with all_w_screen_mirs and add the pigment genes to make fullset
fullset <- full_join(all_w_screen, piggenes, by = "gene_id")
View(fullset)

###Below, between lines of ### is the process I used to get the mature sequences out of the "mirna_alignment"
##I used the output file to clean up the miRNA naming issues in the data. 
## (specifically by searching the sequences of miRNAs to compare to the inserts in my UAS-miRNA flies)
##################################################################################################
##changing mirna_alignment to all lowercase so that it will serve as an identifier of unique
##mirna sequence rather than showing which nucleotides align. that way i can see if any of
##the uniquely named mirnas have identical sequences.
lowercase_miralign <- fullset %>% mutate(mirna_alignment = tolower(mirna_alignment))
head(lowercase_miralign)

uniqetry <- subset(lowercase_miralign, !duplicated(mirna_name))
head(uniqetry)

write.csv(uniqetry, "all_miranda_uniqe_mirnames_only.csv")

mirnames_and_miraligns <- subset(uniqetry, select = c(X.mirbase_acc, mirna_name, mirna_alignment, OE_pigment_PT))
write.csv(mirnames_and_miraligns, "all_miranda_plus_screen_unique_subset.csv")
for_reversal <- read.csv("all_miranda_plus_screen_unique_subset_dashes_removed.csv")
class(for_reversal$mirna_alignment)
for_reversal$mirna_alignment <- as.character(for_reversal$mirna_alignment)
class(for_reversal$mirna_alignment)
head(for_reversal)

strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

with_reversed_aligns <- mutate(for_reversal, mirna_sequence = strReverse(mirna_alignment))

write.csv(with_reversed_aligns, "all_miranda_w_screen_mirna_seuqence.csv")
#################################################################################################

##make lists of unique genes (targets) and mirs (regulators)
gene_list <- unique(fullset$gene_id)
anyNA(gene_list)
##remove NA from gene_list
gene_list <- gene_list[!is.na(gene_list)]
length(gene_list)
mir_list <- unique(fullset$mirna_name)
anyNA(mir_list)
##remove NA from mir_list
mir_list <- mir_list[!is.na(mir_list)]
length(mir_list)


##trimming off some of the excess data from the fullset data for next steps.
trimfullset<- subset(fullset, select = c(mirna_name, gene_id, gene_symbol, Gene, OE_pigment_PT, OE_pigment_PT_fa6, TF, gene_fxn_pigm, kd_pt_fA6))
head(trimfullset)

##get a database of all genes and each unique miRNA that targets it
df <- c(1:200)
df <- cbind(1:200, df)
for (i in gene_list){
  x <- subset(trimfullset, gene_id == i, select = c(gene_id, mirna_name, gene_symbol))
  mirvec <- unique(x$mirna_name)
  entry <- c(toString(i), mirvec)
  length(entry) <- length(df[,1])
  df <- cbind(df, entry)
}
df[,8]
#accidentally made it so each column is named "entry", so I'm changing the colnames to the gene_id
colnames(df) <- df[1,]
#get rid of cols 1&2, which I just used to set up a matrix long enough to accommodate all lengths of new cols 
targdb <- df[-1,-(1:2)]
## in targdb, each column is all the unique mirs that target a particular gene id. colname is gene id
write.csv(targdb, "Mar212018_mirs_by_targ_gene.csv")
#sum(!is.na(targdb[,##column##])) will give the number of seed matches for gene in that column

##seedcounts = count number of unique mirs targeting each gene id
seedcounts <- vector()
for (i in 1:ncol(targdb)){
  seed <- sum(!is.na(targdb[,i]))
  seedcounts <- c(seedcounts, seed)
}
##check that seedcounts is now as long as gene_list, then get a histogram of seedcounts per gene
length(seedcounts) == length(gene_list)
hist(seedcounts)

##get counts of unique miRNAs targeting each gene in an easy to read dataframe
mircountsbygene <- c("gene_id", "gene_symbol","number_of_mirs_targeting", "TF", "gene_fxn_pigm", "kd_pt_fA6")
for(i in gene_list){
  y <- subset(trimfullset, gene_id == i, select = c(gene_id, mirna_name, gene_symbol, TF, gene_fxn_pigm, kd_pt_fA6))
  uniquemirs <- unique(y$mirna_name)
  allmircounts <- length(uniquemirs[!is.na(uniquemirs)])
  mircountsbygene <- rbind(mircountsbygene, c(toString(i), unique(y$gene_symbol), allmircounts, as.character(y$TF[1]), as.character(y$gene_fxn_pigm[1]), as.character(y$kd_pt_fA6[1])))
}
head(mircountsbygene)
tbl_df(mircountsbygene)
##rename columns to the first row names
colnames(mircountsbygene) <- mircountsbygene[1,]
mircountsbygene <- mircountsbygene[-1,]
mircountsbygene <- as.data.frame(mircountsbygene)
mircountsbygene[,3] <- as.numeric(as.character(mircountsbygene[,3]))
write.csv(mircountsbygene, "Mar212018_mircountsbygene.csv")

hist(mircountsbygene$number_of_mirs_targeting, main = "counts of miRNAs predicted to target each gene", 
     xlab = "unique miRNAs predicted to target a gene's 3' UTR", ylab = "number of genes")

##get a database of all miRNAs and their unique target genes
mirtargdb <- c(1:5000)
mirtargdb <- cbind(1:5000, mirtargdb)
for (i in mir_list){
  z <- subset(trimfullset, mirna_name == i, select = c(gene_id, mirna_name, gene_symbol))
  mirvec <- unique(z$gene_id)
  entry <- c(toString(i), mirvec)
  length(entry) <- length(mirtargdb[,1])
  mirtargdb <- cbind(mirtargdb, entry)
}
mirtargdb[,8]
#accidentally made it so each column is named "entry", so I'm changing the colnames to the gene_id
colnames(mirtargdb) <- mirtargdb[1,]
#get rid of columns 1&2, which I just used to set up a matrix long enough to accommodate all lengths of new cols 
mirtargdb <- mirtargdb[-1,-(1:2)]
## in mirtargdb, each column is all the unique mirs that target a particular gene id. colname is gene id
write.csv(mirtargdb, "Mar212018_genes_targeted_by_mir.csv")

##get counts of unique genes targeted by each miRNA in an easy to read dataframe
targcountsbymir <- data.frame(miRNA = character(), number_of_targets = integer(), 
                              OE_phenotype = character(), OE_phenotype_fA6 = character(), stringsAsFactors = FALSE)
for(i in mir_list){
  zz <- subset(trimfullset, mirna_name == i, select = c(mirna_name, gene_id, OE_pigment_PT, OE_pigment_PT_fa6))
  uniquegenes <- unique(zz$gene_id)
  allgenecounts <- length(uniquegenes[!is.na(uniquegenes)])
  targcountsbymir[nrow(targcountsbymir)+1,] <- (c(toString(i), allgenecounts, as.character(zz$OE_pigment_PT[1]), 
                                                  as.character(zz$OE_pigment_PT_fa6[1])))
}
targcountsbymir$number_of_targets <- as.numeric(targcountsbymir$number_of_targets)
head(targcountsbymir)
write.csv(targcountsbymir, "Feb142018_targcountsbymir.csv")
tbl_df(targcountsbymir)

hist(targcountsbymir$number_of_targets, main = "counts of predicted target genes for each miR", 
     xlab = "unique genes predicted to be regulated", ylab = "number of miRNAs")

###testing to see how I can get a list of all genes targeted by each miRNA
mir276a <- subset(trimfullset, mirna_name=='dme-miR-276a')
head(mir276a)
##this will include cases where the same gene is targeted more than once by mir276a
##meaning some genes will show up on this list more than once.

##subset of full dataset only miRNAs sufficient to darken pigmentation
darkmirdb <- filter(trimfullset, OE_pigment_PT == 'B' | OE_pigment_PT == 'D')
tbl_df(darkmirdb)
## "TF" column is where I indicated whether a pigment gene is or isn't a TF with 'y' or 'n'
## There is only an entry in this column for genes from the piggene list
## So summing !is.na(something$TF) will give the number of pigment genes targeted in the list
#sum(!is.na(darkmirdb$TF))

##subset of full dataset only miRNAs sufficient to lighten pigmentation
lightmirdb <- filter(trimfullset, OE_pigment_PT == 'B' | OE_pigment_PT == 'L')
tbl_df(lightmirdb)
#sum(!is.na(lightmirdb$TF))

#sum(!is.na(trimfullset$TF))

##subset of full dataset only miRNAs sufficient to lighten or darken pigmentation
sufficientmirs <- filter(trimfullset, OE_pigment_PT == 'B' | OE_pigment_PT == 'L' | OE_pigment_PT == 'D')
head(sufficientmirs)
sum(!is.na(sufficientmirs$TF))


## for each data subset (sufficientmirs, lightmirdb, darkmirdb, full), what proportion of 
## targeted genes are "pigment genes". Are there more pigment genes targeted by those 
## sufficient to effect pigmentation (in either direction)?
proportion_piggene <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors = FALSE)
proportion_piggene
n <- c('dataset','target_predictions','pigment_gene_targets','proportion')
colnames(proportion_piggene) <- n
proportion_piggene[1,] <- c('full', nrow(trimfullset), sum(!is.na(trimfullset$TF)), 
                            sum(!is.na(trimfullset$TF))/(nrow(trimfullset)))
proportion_piggene[2,] <- c('lighten_or_darken', nrow(sufficientmirs), sum(!is.na(sufficientmirs$TF)), 
                            sum(!is.na(sufficientmirs$TF))/(nrow(sufficientmirs)))
proportion_piggene[3,] <- c('lighten', nrow(lightmirdb), sum(!is.na(lightmirdb$TF)), 
                            sum(!is.na(lightmirdb$TF))/(nrow(lightmirdb)))
proportion_piggene[4,] <- c('darken', nrow(darkmirdb), sum(!is.na(darkmirdb$TF)), 
                            sum(!is.na(darkmirdb$TF))/(nrow(darkmirdb)))
proportion_piggene <- transform(proportion_piggene, target_predictions = as.numeric(target_predictions), 
                                pigment_gene_targets = as.numeric(pigment_gene_targets), proportion = as.numeric(proportion))
proportion_piggene
write.csv(proportion_piggene, "proportion_piggene.csv")
##looks like miRNAs that affect pigmentation are no more likely than any other miRNA to target pigment genes 


####are the proportions of pigment genes targeted by 'sufficient' miRNAs significantly different than
####or significantly higher than those targeted by all miRNAs?
propstat_2side <- function(x)
  prop.test(sum(!is.na(x$TF)),nrow(x), p=sum(!is.na(trimfullset$TF))/(nrow(trimfullset)), alternative = "two.sided")
propstat_greater <- function(x)
  prop.test(sum(!is.na(x$TF)),nrow(x), p=sum(!is.na(trimfullset$TF))/(nrow(trimfullset)), alternative = "greater")
propstat_less <- function(x)
  prop.test(sum(!is.na(x$TF)),nrow(x), p=sum(!is.na(trimfullset$TF))/(nrow(trimfullset)), alternative = "less")

propstat_2side(sufficientmirs) 
## the proportion of pigment genes targeted by miRs sufficient to change pigmentation is significantly different
## from the proportion of pigment genes targeted by all miRs
propstat_greater(sufficientmirs)
## the proportion of pigment genes targeted by miRs sufficient to change pigmentation is not significantly
## greater than the proportion of pigment genes targeted by all miRs
propstat_less(sufficientmirs)
## the proportion of pigment genes targeted by miRs sufficient to change pigmentation is  significantly
## less than the proportion of pigment genes targeted by all miRs

proptest_df <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors = FALSE)
proptesting <- function(x)
  c(deparse(substitute(x)), propstat_2side(x)$p.value, propstat_greater(x)$p.value, propstat_less(x)$p.value)
colnames(proptest_df)<- c("miRNA set", "prop.test 2-sided p-value", "prop.test greater p-value", "prop.test less p-value")
proptest_df[1,]<-proptesting(sufficientmirs)
proptest_df[2,]<-proptesting(darkmirdb)
proptest_df[3,]<-proptesting(lightmirdb)
proptest_df
## proptest_df gives the p values from comparing the proportions of
## [pigmentation genes with seeds targeted by this set of mirs/all genes with seeds targeted by this set of mirs]

## for each data subset (sufficientmirs, lightmirdb, darkmirdb, full), what proportion of 
## targeted genes are "darkening genes". Are there more darkening genes targeted by those 
## sufficient to effect pigmentation (in either direction)?
## lightening mirs (lightmirdb) would logically be expected to target "darkening genes"
dark_expected <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors = FALSE)
dark_expected
m <- c('dataset','target_predictions','darkgenes_targeted','proportion')
colnames(dark_expected) <- m
dark_expected[1,] <- c('full', nrow(trimfullset), length(which(trimfullset$gene_fxn_pigm == 'D' | trimfullset$gene_fxn_pigm == "B")), 
                       length(which(trimfullset$gene_fxn_pigm == 'D' | trimfullset$gene_fxn_pigm == "B"))/(nrow(trimfullset)))

dark_expected[2,] <- c('lighten_or_darken', nrow(sufficientmirs), length(which(sufficientmirs$gene_fxn_pigm == 'D' | sufficientmirs$gene_fxn_pigm == "B")), 
                       length(which(sufficientmirs$gene_fxn_pigm == 'D' | sufficientmirs$gene_fxn_pigm == "B"))/(nrow(sufficientmirs)))
dark_expected[3,] <- c('lighten', nrow(lightmirdb), length(which(lightmirdb$gene_fxn_pigm == 'D' | lightmirdb$gene_fxn_pigm == "B")), 
                       length(which(lightmirdb$gene_fxn_pigm == 'D' | lightmirdb$gene_fxn_pigm == "B"))/(nrow(lightmirdb)))
dark_expected[4,] <- c('darken', nrow(darkmirdb), length(which(darkmirdb$gene_fxn_pigm == 'D' | darkmirdb$gene_fxn_pigm == "B")), 
                       length(which(darkmirdb$gene_fxn_pigm == 'D' | darkmirdb$gene_fxn_pigm == "B"))/(nrow(darkmirdb)))
dark_expected
write.csv(dark_expected, "dark_expected.csv")

## for each data subset (sufficientmirs, lightmirdb, darkmirdb, full), what proportion of 
## targeted genes are "lightening genes". Are there more lightening genes targeted by those 
## sufficient to effect pigmentation (in either direction)?
## darkening mirs (darkmirdb) would logically be expected to target "lightening genes"

light_expected <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors = FALSE)
light_expected
h <- c('dataset','target_predictions','lightgenes_targeted','proportion')
colnames(light_expected) <- h
light_expected[1,] <- c('full', nrow(trimfullset), length(which(trimfullset$gene_fxn_pigm == 'L' | trimfullset$gene_fxn_pigm == "B")), 
                        length(which(trimfullset$gene_fxn_pigm == 'L' | trimfullset$gene_fxn_pigm == "B"))/(nrow(trimfullset)))

light_expected[2,] <- c('lighten_or_darken', nrow(sufficientmirs), length(which(sufficientmirs$gene_fxn_pigm == 'L' | sufficientmirs$gene_fxn_pigm == "B")), 
                        length(which(sufficientmirs$gene_fxn_pigm == 'L' | sufficientmirs$gene_fxn_pigm == "B"))/(nrow(sufficientmirs)))
light_expected[3,] <- c('lighten', nrow(lightmirdb), length(which(lightmirdb$gene_fxn_pigm == 'L' | lightmirdb$gene_fxn_pigm == "B")), 
                        length(which(lightmirdb$gene_fxn_pigm == 'L' | lightmirdb$gene_fxn_pigm == "B"))/(nrow(lightmirdb)))
light_expected[4,] <- c('darken', nrow(darkmirdb), length(which(darkmirdb$gene_fxn_pigm == 'L' | darkmirdb$gene_fxn_pigm == "B")), 
                        length(which(darkmirdb$gene_fxn_pigm == 'L' | darkmirdb$gene_fxn_pigm == "B"))/(nrow(darkmirdb)))
light_expected
write.csv(light_expected, "light_expected.csv")

###################################################################################


##Playing around with adding new info from flybase (expression,etc) but will need to add FBids

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


## select a set of 80 random genes (because there's 80 pigment genes in my list)
## count how many times these genes are regulated by miRNAs that are
## 1) sufficient to darken pigmentation 2) sufficient to lighten pigmentation 3) not sufficient to affect pigmentation
initiate_df <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
initiate_df
count<-1
for (i in 1:100){
  random80 <- sample(miranda_fbgn$fbgn, size = 80, replace = FALSE)
  randoset <- subset(fullset_with_fbgn, fbgn %in% random80 == TRUE, select = c(OE_pigment_PT, Gene))
  pigs <- subset(randoset, !is.na(randoset$Gene), select = Gene)
  pigs <- unique(pigs)
  initiate_df[nrow(initiate_df)+1,]<- c(count,
                                        length(which(randoset$OE_pigment_PT=="D" | randoset$OE_pigment_PT=="B")),
                                        length(which(randoset$OE_pigment_PT=="L" | randoset$OE_pigment_PT=="B")),
                                        length(which(randoset$OE_pigment_PT=="N")),
                                        nrow(pigs))
  count <- count+1
}
boxplot(initiate_df[,2:4])


### for each set of genes (all genes, darkening genes, and lightening genes)
### get a count for the number of seedmatches to darkmirs, lightmirs, and noeffectmirs
pigset <- subset(fullset_with_fbgn, gene_fxn_pigm == "D" | gene_fxn_pigm == "B" | gene_fxn_pigm == "L", select = c(OE_pigment_PT, Gene))
length(!is.na(unique(pigset$Gene)))
head(pigset)
pigsetcount <- data.frame(piggeneset=character(),gene_count=integer(),darkmirs=integer(), lightmirs=integer(), noeffectmirs=integer(), stringsAsFactors = FALSE)
pigsetcount[1,] <- c("all", length(!is.na(unique(pigset$Gene))),
                     length(which(pigset$OE_pigment_PT=="D" | pigset$OE_pigment_PT=="B")),
                     length(which(pigset$OE_pigment_PT=="L" | pigset$OE_pigment_PT=="B")),
                     length(which(pigset$OE_pigment_PT=="N")))


darkset <- subset(fullset_with_fbgn, gene_fxn_pigm == "D" | gene_fxn_pigm == "B", select = c(OE_pigment_PT, Gene))
length(!is.na(unique(darkset$Gene)))
head(darkset)
pigsetcount[2,] <- c("darkening", length(!is.na(unique(darkset$Gene))),
                     length(which(darkset$OE_pigment_PT=="D" | darkset$OE_pigment_PT=="B")),
                     length(which(darkset$OE_pigment_PT=="L" | darkset$OE_pigment_PT=="B")),
                     length(which(darkset$OE_pigment_PT=="N")))

lightset <- subset(fullset_with_fbgn, gene_fxn_pigm == "L" | gene_fxn_pigm == "B", select = c(OE_pigment_PT, Gene))
length(!is.na(unique(lightset$Gene)))
head(lightset)
pigsetcount[3,] <- c("lightening", length(!is.na(unique(lightset$Gene))),
                     length(which(lightset$OE_pigment_PT=="D" | lightset$OE_pigment_PT=="B")),
                     length(which(lightset$OE_pigment_PT=="L" | lightset$OE_pigment_PT=="B")),
                     length(which(lightset$OE_pigment_PT=="N")))
pigsetcount

### I think that the counts of mirs targeting the random set of genes may be lower than 
### the counts of mirs targeting pigment genes because of differences in the composition
### of the gene lists. Specifically, I think there are probably more TFs in the pigggens
### will try to normalize this.

length(which(piggenes$TF == 'y')) ##55
length(which(piggenes$TF == 'n'))  ##25
length(which(darkgenes$TF == 'y')) ##41
length(which(darkgenes$TF == 'n')) ##16
length(which(lightgenes$TF == 'y')) ##26
length(which(lightgenes$TF == 'n')) ##10


### adding column to show which target genes in miranda db are annotated with 
### GO term 0003700 "transcription factor activity, sequence-specific DNA binding"
full_tf_list <- read.delim("GO_0003700_tf_activity_DNA_bind.txt", sep="")
tf <- c("tf")
tf_label <- rep(tf, length.out = nrow(full_tf_list))
full_tf_list[,2]<-tf_label
full_tf_list <- data.frame(full_tf_list, stringsAsFactors = FALSE)
dim(full_tf_list)
head(full_tf_list)
colnames(full_tf_list)<-c("fbgn","is_tf")
fullset_with_fbgn_tf <- left_join(fullset_with_fbgn,full_tf_list,by="fbgn")
length(which(fullset_with_fbgn_tf$is_tf == "tf"))
miranda_fbgn_tf <- left_join(miranda_fbgn,full_tf_list, by="fbgn")
length(which(miranda_fbgn_tf$is_tf == "tf"))

### making sure to have the same number of genes and trying to match proportion of tfs
### permute and select 100 random sets of genes of the same size as the pigmentation gene set
### 3 sets: 
### all pigmentation genes (80 genes, 55tf)
### darkening genes (57 genes, 41tf)
### lightening genes (36 genes, 26tf)

### first up: sample size 80 genes, compare to all pigment genes
permute_tf_condition80 <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf, is_tf =="tf")
  wotf <- filter(miranda_fbgn_tf, is.na(miranda_fbgn_tf$is_tf) == TRUE)
  rando55tf <- sample(wtf$fbgn, size = 55, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 25, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_tf_condition80[nrow(permute_tf_condition80)+1,]<- c(count,
                                                              length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                                              length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                                              length(which(rando$OE_pigment_PT=="N")),
                                                              nrow(pigs))
  count <- count+1
}
permute_tf_condition80
boxplot(permute_tf_condition80[,2:4])

## next sample size 57 genes, compare to darkening genes
permute_tf_condition_dark <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf, is_tf =="tf")
  wotf <- filter(miranda_fbgn_tf, is.na(miranda_fbgn_tf$is_tf) == TRUE)
  rando55tf <- sample(wtf$fbgn, size = 41, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 16, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_tf_condition_dark[nrow(permute_tf_condition_dark)+1,]<- c(count,
                                                                    length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                                                    length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                                                    length(which(rando$OE_pigment_PT=="N")),
                                                                    nrow(pigs))
  count <- count+1
}
permute_tf_condition_dark
boxplot(permute_tf_condition_dark[,2:4])

## next sample size 57 genes, compare to lightening genes
permute_tf_condition_light <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf, is_tf =="tf")
  wotf <- filter(miranda_fbgn_tf, is.na(miranda_fbgn_tf$is_tf) == TRUE)
  rando55tf <- sample(wtf$fbgn, size = 26, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 10, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_tf_condition_light[nrow(permute_tf_condition_light)+1,]<- c(count,
                                                                      length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                                                      length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                                                      length(which(rando$OE_pigment_PT=="N")),
                                                                      nrow(pigs))
  count <- count+1
}
permute_tf_condition_light
boxplot(permute_tf_condition_light[,2:4])

### add genes that are expressed in late pupa or early adult to data
express <- read.delim("FlyBase_expressed_late_pupa_early_adult_moderate.txt", sep="")
on <- c("on")
exp <- rep(on, length.out = nrow(express))
express[,2]<-exp
express <- data.frame(express, stringsAsFactors = FALSE)
dim(express)
head(express)
colnames(express)<-c("fbgn","exp")
fullset_with_fbgn_tf_exp <- left_join(fullset_with_fbgn_tf,express,by="fbgn")
length(which(fullset_with_fbgn_tf_exp$exp == "on"))
miranda_fbgn_tf_exp <- left_join(miranda_fbgn_tf,express, by="fbgn")
length(which(miranda_fbgn_tf_exp$exp == "on"))

############REDO PERMUTATIONS WITH TF AND EXPRESSION
### Now only sampling genes that are expressed in late pupa to early adult AND have
### ratio of TFs comparable to pigment gene set for comparison

### first up: sample size 80 genes, compare to all pigment genes
#####WHAT IS GOING ON WHY ARE THERE ONLY 25 GENES THAT ARE EXPRESSED IN ADULTS AND ARE TF?
permute_80_exp <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf_exp, is_tf =="tf" & exp=="on")
  wotf <- filter(miranda_fbgn_tf_exp, is.na(miranda_fbgn_tf_exp$is_tf) == TRUE & exp=="on")
  rando55tf <- sample(wtf$fbgn, size = 55, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 25, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf_exp, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_80_exp[nrow(permute_80_exp)+1,]<- c(count,
                                              length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                              length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                              length(which(rando$OE_pigment_PT=="N")),
                                              nrow(pigs))
  count <- count+1
}
permute_80_exp
boxplot(permute_80_exp[,2:4])

## next sample size 57 genes, compare to darkening genes
permute_dark_exp <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf_exp, is_tf =="tf" & exp=="on")
  wotf <- filter(miranda_fbgn_tf_exp, is.na(miranda_fbgn_tf_exp$is_tf) == TRUE & exp=="on")
  rando55tf <- sample(wtf$fbgn, size = 41, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 16, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf_exp, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_dark_exp[nrow(permute_dark_exp)+1,]<- c(count,
                                                  length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                                  length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                                  length(which(rando$OE_pigment_PT=="N")),
                                                  nrow(pigs))
  count <- count+1
}
permute_dark_exp
boxplot(permute_dark_exp[,2:4])

## next sample size 57 genes, compare to lightening genes
permute_light_exp <- data.frame(rep=integer(),darkmirs=integer(),lightmirs=integer(),noeffect_mirs=integer(),piggenes=integer())
count <- 1
for (i in 1:100){
  wtf <- filter(miranda_fbgn_tf_exp, is_tf =="tf" & exp=="on")
  wotf <- filter(miranda_fbgn_tf_exp, is.na(miranda_fbgn_tf_exp$is_tf) == TRUE & exp=="on")
  rando55tf <- sample(wtf$fbgn, size = 26, replace = FALSE)
  rando25nottf <- sample(wotf$fbgn, size = 10, replace = FALSE)
  rando80 <- c(rando55tf,rando25nottf)
  rando <- subset(fullset_with_fbgn_tf_exp, fbgn %in% rando80 == TRUE, select = c(OE_pigment_PT, Gene, fbgn))
  pigs <- subset(rando, !is.na(rando$Gene), select = Gene)
  pigs <- unique(pigs)
  permute_light_exp[nrow(permute_light_exp)+1,]<- c(count,
                                                    length(which(rando$OE_pigment_PT=="D" | rando$OE_pigment_PT=="B")),
                                                    length(which(rando$OE_pigment_PT=="L" | rando$OE_pigment_PT=="B")),
                                                    length(which(rando$OE_pigment_PT=="N")),
                                                    nrow(pigs))
  count <- count+1
}
permute_light_exp
boxplot(permute_light_exp[,2:4], main="57 random genes vs light genes")



##### this works - split_IDs$primary_FBgn[split_IDs$gene_symbol == 'y']


