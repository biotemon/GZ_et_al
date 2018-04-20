#Plot community a versitile script to plot stacked bar plots
#Version Nov 17. 2017

(list=ls()) 
graphics.off() 
#Stacked bar for community struture GZ project.

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(plotly)

#Threshold of lowest percentage to show as an individual taxonomy
TH = SETTHRESHOLDHERE

#Set working directory
setwd("SETWORKINGDIRHERE")

#load data
#taxonomyXcounts <- read.delim("/Users/tito-admin/Dropbox/ResearchLisa/Temp14/Methylothrops/methylotrophs_v4.txt")
taxonomyXcounts <- read.delim("SETTAXCOUNTSFILEHERE")
sample_names = c(SETSAMPLENAMESHERE)
samples <- c(1:length(sample_names))
#samples <- c(7, 8, 9, 10, 11, 12) #Only for the 16S data


#simplifyTaxonomyMatrix(my_relative, my_absolute, my_level, my_taxonomycounts)
#my_level ==> Kingdom level 1
#my_level ==> Phylum level 2
#my_level ==> Class level 3
#my_level ==> Order level 4
#my_level ==> Family level 5
#my_level ==> Genus level 6
#my_level ==> Species level 7

#simplifyTaxonomyMatrix(relative_abundance2, tax_counts2, 7, taxonomyXcounts)

#Find unique headers
taxones = c()
for (i in samples)
{
  taxonomyXcounts_sample <- filter(taxonomyXcounts, ASSEMBLY_ID == i)
  sums_by_taxon <- ddply(taxonomyXcounts_sample, .(TAX_ID), summarise, READ_COUNTS=sum(READ_COUNTS))
  vec <- as.vector(sums_by_taxon$TAX_ID)
  taxones = c(taxones, vec)
}

uniq_taxons = unique(taxones)

# Initialize matrix

n = length(samples)
m = length(uniq_taxons)

tax_counts = matrix(rep(0,n*m),nrow=n,ncol=m)


#Fill the matrix
for (i in 1:length(samples))
{
  taxonomyXcounts_sample <- filter(taxonomyXcounts, ASSEMBLY_ID == samples[i])
  sums_by_taxon <- ddply(taxonomyXcounts_sample, .(TAX_ID), summarise, READ_COUNTS=sum(READ_COUNTS))
  vec <- as.vector(sums_by_taxon$TAX_ID)
  vec2 <- as.vector(sums_by_taxon$READ_COUNTS)
  uniq_taxons[which(is.na(uniq_taxons))] <- 0
  vec[which(is.na(vec))] <- 0
  for (j in 1:length(vec)){
    for(k in 1:length(uniq_taxons)){
      if (vec[j] == uniq_taxons[k]){
        #outstring  = paste(vec[j], uniq_taxons[k], vec2[j], sep = " ")
        #print(outstring) 
        tax_counts[i,k] <- as.numeric(vec2[j])
      }
    }
    
  }
}

#if(uniq_taxons[1] == ""){
#  uniq_taxons[1] = "Others"
#}

#Relative Abundance
SumSamples = rowSums(tax_counts)
relative_abundance = (tax_counts *100) / SumSamples

colnames(tax_counts) <- uniq_taxons
colnames(relative_abundance) <- uniq_taxons

#RUN HERE THE FUNCTION 

#Functions

#simplifyTaxonomyMatrix <- function(my_relative, my_absolute, my_level, my_taxonomycounts)
#{

my_relative = relative_abundance
my_absolute = tax_counts
my_level = 8
my_taxonomycounts = taxonomyXcounts


#We want to coalesce or merge those taxonomies with lower than a certain threshold (0.5%) into a higher rank of taxonomy.

#Subset higher ranks
high_ranks = setNames(data.frame(my_taxonomycounts[,2], my_taxonomycounts[,7], my_taxonomycounts[,8], my_taxonomycounts[,9], my_taxonomycounts[,10], my_taxonomycounts[,11], my_taxonomycounts[,12], my_taxonomycounts[,13], my_taxonomycounts[,14]), c("ID", "SUPERKINGDOM","KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"))
#remove repeated rows
high_ranks = high_ranks[!duplicated(high_ranks), ]
high_ranks = high_ranks[order(high_ranks$SUPERKINGDOM, high_ranks$KINGDOM, high_ranks$PHYLUM, high_ranks$CLASS, high_ranks$ORDER, high_ranks$FAMILY, high_ranks$GENUS), ]
#Set 0 in NA column as well as for empty taxons
if(any(is.na(high_ranks$ID))){
high_ranks[which(is.na(high_ranks$ID)),1] <- 0
}

#Extend last word to replace NAs and '' cells
for(i in 1:dim(high_ranks)[1]){
  last_term = ''
  for(j in 2:dim(high_ranks)[2]){
    if(!is.na(high_ranks[i,j]))
      {
      if(high_ranks[i,j] != ""){
        if( !grepl("[[:space:]]+", high_ranks[i,j])){
          last_term <- as.character(high_ranks[i,j])
        }
      }
    }
   if(last_term != ''){
     #print(last_term)
      my_levels <- levels(high_ranks[,j])
      my_levels[length(my_levels) + 1 ] <- last_term
      my_levels <- unique(my_levels)
      high_ranks[,j] <- factor(high_ranks[,j], levels = my_levels)
      high_ranks[i,j] <- last_term
   }
    
  }
} 


zerocols <- high_ranks[which(high_ranks[,2] == '' & high_ranks[,3] == '' & high_ranks[,4] == '' & high_ranks[,5] == '' & high_ranks[,6] == '' & high_ranks[,7] == '' & high_ranks[,8] == '' & high_ranks[,1] != 0),1]
zerovec <- as.vector(which(high_ranks[,2] == '' & high_ranks[,3] == '' & high_ranks[,4] == '' & high_ranks[,5] == '' & high_ranks[,6] == '' & high_ranks[,7] == '' & high_ranks[,8] == '' & high_ranks[,1] != 0))
if(any(zerovec)){
high_ranks <- high_ranks[-zerovec,]
}
if(any(zerocols)){
for(i in zerocols){
  colnames(my_absolute)[colnames(my_absolute) == i] <- 0
  colnames(my_relative)[colnames(my_relative) == i] <- 0
}
}





#Loop to coalesce columns names
my_level = my_level - 1
for(k in my_level:2)
{
#k = 2
  matrix_col_names = as.vector(colnames(my_relative))
  
  for(i in 1:dim(my_relative)[2])
  {
    #i = 8
    low_counter = 0
    for(j in 1:dim(my_relative)[1])
    {
      if(my_relative[j,i] < TH)
      {
        low_counter = low_counter + 1
      } 
    }
    if(low_counter == 6)
    {
      #Get the row in taxonomy where is the ID in the ID column.
      if(any(grep("^[0-9]+$",matrix_col_names[i], perl = TRUE, value=FALSE))){
        row_val = which(high_ranks[,1] == matrix_col_names[i])
        #Look for the 
        higher_name = as.character(high_ranks[row_val[1],k])
      }else{
        #higher_name = as.character(filter(high_ranks, CLASS == matrix_col_names[i])[[2]])
        my_k = k+1
        row_val = which(high_ranks[,my_k] == matrix_col_names[i])
        #Look for the 
        higher_name = as.character(high_ranks[row_val[1],k])
      }
      #Change the column ID
      if(!is.na(higher_name)){
        if(higher_name != ''){
          colnames(my_relative)[colnames(my_relative) == matrix_col_names[i]] <- higher_name
          colnames(my_absolute)[colnames(my_absolute) == matrix_col_names[i]] <- higher_name
        }
      }
      
    }
 }
  
  #Merging columns by the same name
  
  my_absolute = t(rowsum(t(my_absolute), colnames(my_absolute)))
  my_relative = t(rowsum(t(my_relative), colnames(my_relative)))
  
}

#Change names of genus
matrix_col_names = as.vector(colnames(my_relative))
for(i in 1:dim(my_relative)[2])
{
  if(any(grep("^[0-9]+$",matrix_col_names[i], perl = TRUE, value=FALSE))){
    row_val = which(high_ranks[,1] == matrix_col_names[i])
    higher_name = as.character(high_ranks[row_val[1],7])
    part2_name = as.character(high_ranks[row_val[1],8])
    if(!is.na(part2_name)){
    if(part2_name != ""){
    if(part2_name != higher_name){
      higher_name = paste(higher_name, part2_name, sep=' ')
    }}}
    colnames(my_relative)[colnames(my_relative) == matrix_col_names[i]] <- higher_name
    colnames(my_absolute)[colnames(my_absolute) == matrix_col_names[i]] <- higher_name
    
  }
}


#Merging columns by the same name... one more time just in case

my_absolute = t(rowsum(t(my_absolute), colnames(my_absolute)))
my_relative = t(rowsum(t(my_relative), colnames(my_relative)))

#}


###### Preparing Matrices for plotting ########

simple_absolute_matrix_2 = my_absolute
simple_relative_matrix_2 = my_relative

#Prepare Matrices for ploting
#Next line binds + create dataframe + keep numeric as numeric and dont change to text
simple_absolute_matrix_2 <- cbind.data.frame(sample_names, simple_absolute_matrix_2)
simple_relative_matrix_2 <- cbind.data.frame(sample_names, simple_relative_matrix_2)

#Reorder Columns by taxonomy

guide_vector = c()
for(i in 1:dim(high_ranks)[1])
  {
  for(j in 2:7)
    {
    my_value = as.character(high_ranks[i,j])
      if(is.na(match(my_value, guide_vector)) && (my_value != ""))
        {
          guide_vector = c(guide_vector, my_value)
        }
    }
  }

#Ordering the tables
headers_in_simple = colnames(simple_absolute_matrix_2)

#Remove 'sample_names' and 'V1' from headers_in_simple

headers_in_simple = headers_in_simple[headers_in_simple != "sample_names"]
headers_in_simple = headers_in_simple[headers_in_simple != "V1"]

#First we order the names later we append columns 
numeric_index=c()
names_index=c()

for(i in 1:length(headers_in_simple)){
  x_name <- as.vector(strsplit(headers_in_simple[i], " "))[[1]][1]
  if(any(grep(x_name, guide_vector, perl = TRUE, value=FALSE))){
    x_val <- as.numeric(grep(x_name, guide_vector, perl = TRUE, value=FALSE))
    numeric_index = c(numeric_index, x_val[1])
    names_index = c(names_index, headers_in_simple[i])
  }else{
    x_val = 1000000
    numeric_index = c(numeric_index, x_val)
    names_index = c(names_index, headers_in_simple[i])
  }
}

the_index_df = cbind.data.frame(numeric_index, names_index)
the_index_df = the_index_df[order(numeric_index),]


#READ THE NEXT THREE LINES#

suited_guide_vector <- as.vector(the_index_df$names_index)
if(colnames(simple_absolute_matrix_2)[2] == "V1"){
suited_guide_vector <- as.vector(c(c("sample_names", "V1"), suited_guide_vector)) #Run this if there is a NA / 0 column in the high ranks
}else{
  suited_guide_vector <- as.vector(c(c("sample_names"), suited_guide_vector))
}



simple_absolute_matrix_3 = simple_absolute_matrix_2[suited_guide_vector]
simple_relative_matrix_3 = simple_relative_matrix_2[suited_guide_vector]

colnames(simple_absolute_matrix_3)[colnames(simple_absolute_matrix_3) == "V1"] <- "Others"
colnames(simple_relative_matrix_3)[colnames(simple_relative_matrix_3) == "V1"] <- "Others"

#Now dataframes are simplified and ordered
#Reshaping data frame
#Making the dataframe from wide to log format

simple_absolute_melt <- data.table::melt(simple_absolute_matrix_3, id.vars='sample_names')
simple_relative_melt <- data.table::melt(simple_relative_matrix_3, id.vars='sample_names')

#Change some names
#Rhodospirillaceae
#AEGEAN-169_marine_group
old_terms_vec <- c("Candidatus_Pelagibacter ubique", 
                   "Candidatus_Pelagibacter ubique HIMB083", 
                   "SAR116_clade",
                   "SAR116_clade SAR116 cluster alpha proteobacterium HIMB100", 
                   "Candidatus_Puniceispirillum marinum", 
                   "SAR11_clade",
                   "Deep_1",
                   "Surface_1",
                   "OM43_clade", 
                   "SAR406_clade", 
                   "Candidatus_Pelagibacter ubique HTCC1002", 
                   "Rhodobacterales Rhodobacterales bacterium HTCC2255",
                   "Rhodobacterales Rhodobacterales bacterium Y4I", 
                   "Rhodobacteraceae Rhodobacteraceae bacterium HTCC2083", 
                   "Rhodobacteraceae Rhodobacteraceae bacterium HTCC2150", 
                   "Rhodobacteraceae Rhodobacteraceae bacterium KLH11", 
                   "Candidatus_Puniceispirillum", 
                   "Candidatus Puniceispirillum marinum IMCC1322",
                   "SAR11_clade uncultured SAR11 cluster alpha proteobacterium H17925_38M03", 
                   "SAR11_clade uncultured SAR11 cluster alpha proteobacterium H17925_48B19", 
                   "SAR11_clade uncultured SAR11 cluster bacterium HF0010_09O16", 
                   "SAR11_clade uncultured SAR11 cluster bacterium HF0770_37D02", 
                   "SAR11_clade uncultured SAR11 cluster bacterium HF4000_37C10", 
                   "Chesapeake−Delaware_Bay", 
                   "LD12_freshwater_group", 
                   "Surface_2", 
                   "Surface_3", 
                   "Surface_4",
                   "SAR11_clade SAR11 cluster bacterium PRT-SC02")

new_terms_vec <- c("Candidatus Pelagibacter ubique", 
                   "Candidatus Pelagibacter ubique HIMB083", 
                   "SAR116 clade", 
                   "SAR116 clade: str. HIMB100", 
                   "Candidatus Puniceispirillum marinum", 
                   "SAR11 clade", 
                   "SAR11 clade: Deep_1", 
                   "SAR11 clade: Surface_1", 
                   "Methylophilales: OM43 clade", 
                   "SAR406 clade", 
                   "Candidatus Pelagibacter ubique HTCC1002", 
                   "Rhodobacterales: str. HTCC2255", 
                   "Rhodobacterales: str. Y4I", 
                   "Rhodobacteraceae: str. HTCC2083", 
                   "Rhodobacteraceae: str. HTCC2150", 
                   "Rhodobacteraceae: str. KLH11", 
                   "Candidatus Puniceispirillum", 
                   "Candidatus Puniceispirillum marinum IMCC1322", 
                   "SAR11 clade: str. H17925_38M03", 
                   "SAR11 clade: str. H17925_48B19", 
                   "SAR11 clade: str. HF0010_09O16", 
                   "SAR11 clade: HF0770_37D02", 
                   "SAR11 clade: HF4000_37C10", 
                   "SAR11 clade: Chesapeake−Delaware Bay", 
                   "SAR11 clade: LD12 freshwater group", 
                   "SAR11 clade: Surface_2", 
                   "SAR11 clade: Surface_3", 
                   "SAR11 clade: Surface_4",
                   "SAR11 clade: str. PRT-SC02")

for(i in 1:length(old_terms_vec)){
  old_term <- old_terms_vec[i]
  term_to_change <- new_terms_vec[i]
  change_val = which(simple_absolute_melt[,2] == old_term)
  #print(paste(old_term,term_to_change,  change_val ))
  
  my_levels <- levels(simple_absolute_melt[,2])
  my_levels[length(my_levels) + 1] <- term_to_change
  my_levels <- unique(my_levels)
  simple_absolute_melt[, 2] <- factor(simple_absolute_melt[, 2], levels = my_levels)
  simple_relative_melt[, 2] <- factor(simple_relative_melt[, 2], levels = my_levels)
  simple_absolute_melt[change_val, 2] <- term_to_change
  simple_relative_melt[change_val, 2] <- term_to_change
  my_levels[which(my_levels == old_term)] <- term_to_change
  my_levels <- unique(my_levels)
  simple_absolute_melt[, 2] <- factor(simple_absolute_melt[, 2], levels = my_levels)
  simple_relative_melt[, 2] <- factor(simple_relative_melt[, 2], levels = my_levels)
}


#Colors
xval =  dim(simple_relative_matrix_3)[2] - 1
#colfunc <- colorRampPalette(brewer.pal(11,"Spectral"))
colfunc <- colorRampPalette(c("#c7c5c4", "#9001c9", "#cc12c0", "#9E0142", "#ff9e9e", "#ff2900", "#ff6f00", "#ffaa00" , "#ffe100",  "#eadf60", "#8db500",  "#b5d153", "#118200", "#72996c", "#00ffd8", "#00b2ff", "#97c6bf",  "#43a898", "#357089",  "#1b44d6", "#849fff", "#5a6384" ))
simple_color_vec <- colfunc(xval)
#simple_color_vec[] <- "#c7c5c4" #Grey
#simple_color_vec[1] <- "#9001c9" #Violet
#simple_color_vec[2] <- "#cc12c0" #Fucsia
#simple_color_vec[3] <- "#9E0142" #Red
#simple_color_vec[5] <- "#ff9e9e" #Light red
#simple_color_vec[6] <- "#ff2900" #coral red
#simple_color_vec[7] <- "#ff6f00" #Orange
#simple_color_vec[8] <- "#ffaa00" #Orange
#simple_color_vec[9] <- "#ffe100" #Yellow
#simple_color_vec[10] <- "#eadf60" #light Yellow
#simple_color_vec[11] <- "#8db500" #Verde mango biche 
#simple_color_vec[12] <- "#b5d153" #Verde menta pastel
#simple_color_vec[13] <- "#118200" #Verde pasto
#simple_color_vec[14] <- "#72996c" #Verde gris
#simple_color_vec[15] <- "#00ffd8" #light blue
#simple_color_vec[16] <- "#00b2ff" #blue
#simple_color_vec[17] <- "#97c6bf" #Light grey blue
#simple_color_vec[18] <- "#43a898" #grey blue
#simple_color_vec[19] <- "#357089" #dark grey blue
#simple_color_vec[20] <- "#1b44d6" #royal blue
#simple_color_vec[21] <- "#849fff" #royal light blue
#simple_color_vec[22] <- "#5a6384" #royal grye blue

pdf("GZ_MGRAST_abs_cutoff_2.pdf", width=12, height=7)
ggplot(data = simple_absolute_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity") + theme_classic() + theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black")) + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1))
dev.off()

pdf("GZ_MGRAST_rel_cutoff_2.pdf", width=10, height=8)
ggplot(data = simple_relative_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(colour="black", stat = "identity") + theme_classic() + theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black")) +  scale_colour_manual("black") + scale_fill_manual(values = simple_color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=1))
dev.off()


















############## PLOTS WITH THE WHOLE TAXON LIST ##########

#Reshaping data frame
#Making the dataframe from wide to log format

#long_DF <- tax_counts2 %>% gather(Taxonomy, ReadCounts, 2:(m+1))
tax_counts_melt <- data.table::melt(tax_counts2, id.vars='sample_names')

#Preparing palette
colfunc <- colorRampPalette(brewer.pal(11,"Spectral"))
color_vec <- colfunc(m)
my_color_vec = c("#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                    "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#AD1245", "#ff7c11", "#9E0142",                  "#A50943", "#AD1245", "#ffe311", "#9E0142", "#A50943",                  "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#AD1245", "#B51A47", "#BC2349",                  "#9E0142", "#A50943", "#a5ed63", "#55B0AC", "#4EA7AF",                  "#388FB9", "#3287BC", "#387FB8", "#3F77B5", "#456FB1",                 "#4B67AD", "#515FA9", "#5757A5", "#5E4FA2", "#456FB1",                 "#388FB9", "#3287BC", "#387FB8", "#3F77B5", "#456FB1",                 "#4B67AD", "#515FA9", "#5757A5", "#5E4FA2", "#456FB1",                 "#4B67AD", "#515FA9")

#Using ggplot
pdf("GZ_CommuStruc_Phylum_2.pdf", width=12, height=9)
ggplot(data = tax_counts_melt, aes(x = sample_names, y = value, fill = variable)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black")) + scale_fill_manual(values = color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=3))
dev.off()

#using plotly
p <-  plot_ly(tax_counts_melt, x= ~sample_names, y = ~value, type = 'bar', name = ~variable, color = ~variable) %>% layout (yaxis = list(title = 'Counts'), barmode = 'stack')




#Reshaping data frame
#Making the dataframe from wide to log format

tax_counts_melt_rel <- data.table::melt(relative_abundance2, id.vars='sample_names')

pdf("GZ_CommuStruc_Phylum_rel_2.pdf", width=12, height=9)
ggplot(data = tax_counts_melt_rel, aes(x = sample_names, y = value, fill = variable)) + geom_bar(stat = "identity") + theme_classic() + theme(axis.text.x = element_text(color="black"),axis.text.y = element_text(color="black")) + scale_fill_manual(values = color_vec) + scale_y_continuous(name="Read Counts", labels = scales::comma, expand = c(0, 0)) + guides(fill=guide_legend(ncol=3))
dev.off()


#################### ---- 


#Generating simplified plots

simple_relative_matrix = relative_abundance2
simple_absolute_matrix = tax_counts2

#We want to coalesce or merge those taxonomies with lower than a certain threshold (0.5%) into a higher rank of taxonomy.

#Subset higher ranks

high_ranks = setNames(data.frame(taxonomyXcounts$KINGDOM, taxonomyXcounts$PHYLUM, taxonomyXcounts$CLASS), c("KINGDOM", "PHYLUM", "CLASS"))
#remove repeated rows
high_ranks = high_ranks[!duplicated(high_ranks), ]
high_ranks = high_ranks[order(high_ranks$KINGDOM, high_ranks$PHYLUM), ]
#Check that the next vector has the form "sample_name" "Others" ...
matrix_col_names = as.vector(colnames(relative_abundance2))

for(i in 3:dim(relative_abundance2)[2]){
  low_counter = 0
  for(j in 1:dim(relative_abundance2)[1]){
   if(relative_abundance2[j,i] < TH){
     low_counter = low_counter + 1
   } 
  }
  if(low_counter == 6){
    
    higher_name = as.character(filter(high_ranks, CLASS == matrix_col_names[i])[[2]])
    #Change the column name
    names(simple_relative_matrix)[names(simple_relative_matrix) == matrix_col_names[i]] <- higher_name
    names(simple_absolute_matrix)[names(simple_absolute_matrix) == matrix_col_names[i]] <- higher_name
  }
}

#Before adding columns we need to remove the first column which is text
simple_absolute_matrix$sample_names <- NULL
simple_relative_matrix$sample_names <- NULL

#Merging columns by the same name

simple_absolute_matrix_2 = t(rowsum(t(simple_absolute_matrix), colnames(simple_absolute_matrix)))
simple_relative_matrix_2 = t(rowsum(t(simple_relative_matrix), colnames(simple_relative_matrix)))

#Prepare Matrices for ploting
#Next line binds + create dataframe + keep numeric as numeric and dont change to text
simple_absolute_matrix_3 <- cbind.data.frame(sample_names, simple_absolute_matrix_2)
simple_relative_matrix_3 <- cbind.data.frame(sample_names, simple_relative_matrix_2)

#Reorder Columns by taxonomy
#First get a guide vector
high_ranks_v2 = filter(high_ranks, KINGDOM != "" & PHYLUM != "" & CLASS != "")

#Remove duplicated lines
high_ranks_v3 <- high_ranks_v2[-which(as.character(high_ranks_v2$KINGDOM) == as.character(high_ranks_v2$PHYLUM)), ]
high_ranks_v4 <- high_ranks_v3[-which(as.character(high_ranks_v3$PHYLUM) == as.character(high_ranks_v3$CLASS)), ]
rm(high_ranks_v2,high_ranks_v3)

guide_vector = c(as.character(high_ranks_v4[1,1]))
for(i in 1:dim(high_ranks_v4)[1]){
  for(j in 1:dim(high_ranks_v4)[2]){
    my_value = as.character(high_ranks_v4[i,j])
    if(j == 1){
      if(is.na(match(my_value, guide_vector))){
        guide_vector = c(guide_vector, my_value)
      }
    }
    else if(j == 2){
      if(is.na(match(my_value, guide_vector))){
        guide_vector = c(guide_vector, my_value)
      }
    }
    else{
      guide_vector = c(guide_vector, my_value)
    }
    
  }
}

#Ordering the tables
headers_in_simple = colnames(simple_absolute_matrix_3)
simple_absolute_matrix_4 = as.data.frame(simple_absolute_matrix_3$sample_names)
simple_relative_matrix_4 = as.data.frame(simple_relative_matrix_3$sample_names)
names_4 = c("sample_names")
#This for append columsn in order
for(i in 1:length(guide_vector)){
  if(!is.na(match(guide_vector[i], headers_in_simple))){
    simple_absolute_matrix_4 = cbind.data.frame(simple_absolute_matrix_4, simple_absolute_matrix_3[,guide_vector[i]])
    simple_relative_matrix_4 = cbind.data.frame(simple_relative_matrix_4, simple_relative_matrix_3[,guide_vector[i]])
    names_4 = c(names_4, guide_vector[i]) 
  }
}
colnames(simple_absolute_matrix_4) <- names_4
colnames(simple_relative_matrix_4) <- names_4

