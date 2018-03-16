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
