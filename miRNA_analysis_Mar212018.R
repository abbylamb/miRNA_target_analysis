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


