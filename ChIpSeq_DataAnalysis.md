Class\_exercise
================
Christopher Ebuka Ojukwu
3/17/2023


## First I will read in each replicate file

``` r
#link to my path
broadpeakfilepath <- "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/group_folder/group/results/bwa/mergedLibrary/macs/broadPeak"

#imports my peaks into a list
peak_list <- import_peaks(consensus_file_path = broadpeakfilepath)

#this counts the number of peaks for each individual replicate
peak_num <- sapply(peak_list, length) %>% as.data.frame()

#This labels the column for the number of peaks
names(peak_num) <- c("num_peaks")

#this names the columns for dbp and replicate
peak_num <-peak_num %>%
  rownames_to_column(var = "dbp") %>%
  separate(col = dbp,  into = c('dbp', 'replicate'), sep = "_")

peak_num
```

    ##      dbp replicate num_peaks
    ## 1  CEBPZ        R1       273
    ## 2  CEBPZ        R2       402
    ## 3   CHD2        R1     22119
    ## 4   CHD2        R2     13012
    ## 5   CTCF        R1     36870
    ## 6   CTCF        R2     65526
    ## 7   CTCF        R3     83141
    ## 8   CTCF        R4     55484
    ## 9   CTCF        R5     52373
    ## 10  CTCF        R6     44241
    ## 11  CTCF        R7     71705
    ## 12  CTCF        R8     68383
    ## 13  ELF1        R1     25895
    ## 14  ELF1        R2     35371
    ## 15  ELF1        R3     23188
    ## 16  ELF1        R4     15436
    ## 17 EP300        R1     38462
    ## 18 EP300        R2     59348
    ## 19 EP300        R3      4471
    ## 20 EP300        R4      6888
    ## 21 EP300        R5     10270
    ## 22 EP300        R6       882

##Now I am going to create consensus peaks for each protein

``` r
#this gets the name of dbps (I have 5 of them).
dbps <- unique(sapply(names(peak_list), function(x){
   unlist(strsplit(x, "_"))[1]
}))

#this finds all the unique peaks in all the replicates
consensus_list <- lapply(dbps, consensus_from_reduced, peak_list)
names(consensus_list) <- dbps

#Adding names to dbps
names(consensus_list) <- dbps 

# Making a data frame of the consensus peaks for each DBP
num_consensus_peaks <- sapply(consensus_list, length) %>% 
  as.data.frame() %>%
  rownames_to_column( var = "dbp") %>%
  dplyr::rename(number_consensus_peaks = ".")

# Merging the consensus peaks into the peak_numbers data frame
peak_num <- left_join(peak_num, num_consensus_peaks)
```

``` r
# Printing table
peak_num
```

    ##      dbp replicate num_peaks number_consensus_peaks
    ## 1  CEBPZ        R1       273                    172
    ## 2  CEBPZ        R2       402                    172
    ## 3   CHD2        R1     22119                   8815
    ## 4   CHD2        R2     13012                   8815
    ## 5   CTCF        R1     36870                  28799
    ## 6   CTCF        R2     65526                  28799
    ## 7   CTCF        R3     83141                  28799
    ## 8   CTCF        R4     55484                  28799
    ## 9   CTCF        R5     52373                  28799
    ## 10  CTCF        R6     44241                  28799
    ## 11  CTCF        R7     71705                  28799
    ## 12  CTCF        R8     68383                  28799
    ## 13  ELF1        R1     25895                  11874
    ## 14  ELF1        R2     35371                  11874
    ## 15  ELF1        R3     23188                  11874
    ## 16  ELF1        R4     15436                  11874
    ## 17 EP300        R1     38462                    568
    ## 18 EP300        R2     59348                    568
    ## 19 EP300        R3      4471                    568
    ## 20 EP300        R4      6888                    568
    ## 21 EP300        R5     10270                    568
    ## 22 EP300        R6       882                    568

``` r
# the path for exporting my consensus peaks
exportpath <- "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_consensus_peaks/"

#exports all the unique peaks in all the replicates
# export consensus peaks to results folder
for(i in 1:length(consensus_list)) {
rtracklayer::export(consensus_list[[i]], paste0(exportpath, names(consensus_list)[i], "_consensus_peaks.bed") )}
```

## Now I am going to make my consensus peaks compatible with UCSC genome browser

``` r
#makes a file list
consensus_file_list <- list.files("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_consensus_peaks/", full.names = T, pattern = ".bed")

#adding column names
peaks <- lapply(consensus_file_list, read.table, col.names = c("chr", "start", "end", "name", "score", "strand"))
names(peaks) <- dbps
# double check order by looking at consensus_file_list is same order as dbps

# make chromosomes of interest object
canonical_chr <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY")

# let's use lapply with filter funciton to cannonical_chr
peaks <- lapply(peaks, function(x) x %>% filter(chr %in% canonical_chr))

# export to consensus peaks folder
new_filenames <- paste0("my_results/my_consensus/", names(peaks), "_consensus.bed")

for(i in 1:length(peaks)) {
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}


# we are using paste0 to print the header text
# then we are adding the names as the value the header equals with 'names' function
headers <- paste0("track type=bed name=", names(peaks))

# print out consensus peak files in a results/UCSC directory
# creating a path to export after we add header in for loop below
new_filenames <- paste0("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_ucsc_consensus/", names(peaks), ".bed")

for(i in 1:length(peaks)) {
  #Write the header line
  writeLines(headers[[i]], new_filenames[[i]])
  #Append the broadPeak table data
  
  write.table(peaks[[i]], new_filenames[[i]],
              sep = "\t", col.names = FALSE, row.names = FALSE,
              quote = FALSE, append = TRUE)
}  
```

## I am curious if my proteins are transcription factors so I will use the annotations
## in a cell paper I found and see

``` r
# import consensus bed files
consensusPeakPath <- "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_consensus/"

consensus_peaks_files <- list.files(consensusPeakPath, 
                                             pattern = "*consensus.bed",
                                             full.names = TRUE)

# lapply with import function to make a list of GRanges
consensus_peaks <- lapply(consensus_peaks_files, rtracklayer::import)

names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_consensus//|_consensus.bed","", consensus_peaks_files)

# import gencode Granges
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf")

# loading in the number of peaks each DBP has -- using length.
num_peaks_df <- data.frame("dbp" = names(consensus_peaks),
                           "num_peaks" = sapply(consensus_peaks, length))

# total amount of the genome covered by all the peaks for a given DBP.
num_peaks_df$total_peak_length <- sapply(consensus_peaks, function(x) sum(width(x)))


url <- "https://www.cell.com/cms/10.1016/j.cell.2018.01.029/attachment/ede37821-fd6f-41b7-9a0e-9d5410855ae6/mmc2.xlsx"

destination_for_url <- "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/TF_annotations.xlsx"

# to download we can use download.file
download.file(url, destination_for_url)


#redx1::read_excel to import
human_tfs <- readxl::read_excel("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/TF_annotations.xlsx",
                                sheet = 2, skip = 1)
```

``` r
# let's rename the 4th column to indicate if it is a TF.
names(human_tfs)[4] <- "is_tf"

#none of my dbps have entries in here
length(which(tolower(num_peaks_df$dbp) %in% tolower(human_tfs$Name)))
```

    ## [1] 5

``` r
# first let's filter and grab the first 4 columns that match DBPs in num_peaks_df
human_tfs <- human_tfs[tolower(human_tfs$Name) %in% tolower(num_peaks_df$dbp), 1:4]


# adding new column names
names(human_tfs) <- c("ensembl_id",
                      "dbp",
                      "dbd",
                      "tf")

# Now comes the actual merge. Let's look at it first
# ?merge

num_peaks_df <- merge(num_peaks_df, human_tfs, all.x = T)

# Let's check how many NAs -- we should have some missing values.
dim(num_peaks_df[is.na(num_peaks_df$tf),])
```

    ## [1] 0 6

``` r
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/num_peaks_df.csv")

num_peaks_df
```

    ##     dbp num_peaks total_peak_length      ensembl_id      dbd  tf
    ## 1 CEBPZ      1892            740707 ENSG00000115816  Unknown Yes
    ## 2  CHD2     96965          49710914 ENSG00000173575 Myb/SANT  No
    ## 3  CTCF    316789         204671555 ENSG00000102974  C2H2 ZF Yes
    ## 4  ELF1    130614          83297027 ENSG00000120690      Ets Yes
    ## 5 EP300      6248           5503663 ENSG00000100393  Unknown  No

CEBPZ, CTCF, and ELF1 are all transcriptions factors.

## Now I want to compare a protein with a previous analysis

``` r
# goto UCSC genome browser and load in a peak file for a given protein
# load in the data for the same protein from the previous analysis
# compare how your consensus peaks are similar or different to previous analyses

knitr::include_graphics("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/ucsc_cebpz_shot.jpg")
```

<img src="christopher_new_files/figure-gfm/ucsc_cebpz_shot.jpg" width="1888" />

``` r
#The peaks from my protein, CEBPZ, are similar to the analysis from previous years.
```

## Now I am going to determine how my peaks for each protein overlap annotations of the genome

## First I will find the overlaps between my consensus peaks with promoters of lncRNA and mRNA promoters

``` r
#we could get just the genes
# now just type == gene annotations
gencode_genes <- gencode_gr[gencode_gr$type == "gene"] 
table(gencode_gr$type)
```

    ## 
    ##           gene     transcript           exon            CDS    start_codon 
    ##          60609         227462        1372308         761508          87662 
    ##     stop_codon            UTR Selenocysteine 
    ##          79913         310193            119

``` r
# exporting just the genes file (we will save all the .Rdata too at the end)
#also saves time; instead of reloading gencode_gr everytime
rtracklayer::export(gencode_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/gencode_genes.gtf")

#maybe redundant
# mRNA genes (called "protein_coding") in this version of gencode changes sometimes !
mrna_genes <- gencode_genes[gencode_genes$gene_type %in% "protein_coding"] 

rtracklayer::export(mrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/mrna_genes.gtf")

# now doing a second index for lncRNA:
lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% "lncRNA"] 

rtracklayer::export(lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/lncrna_genes.gtf")

# both mRNA and lncRNA annotations together.
mrna_lncrna_genes <- gencode_genes[gencode_genes$gene_type %in% c("protein_coding","lncRNA")]
rtracklayer::export(mrna_lncrna_genes, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/mrna_lncrna_genes.gtf")

#we can now import this in the future
# starting annotation file that we will use moving forward.
lncrna_mrna_genes <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/mrna_lncrna_genes.gtf")

# Nice that loaded so much faster -- lets see as DF
lncrna_mrna_genes_df <- lncrna_mrna_genes %>% as.data.frame()

# creating and exporting promoter annotations
lncrna_mrna_promoters <- promoters(lncrna_mrna_genes, upstream = 1000, downstream = 1000)

#feature for future read in
rtracklayer::export(lncrna_mrna_promoters, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/my_results/my_gene_annotation/lncrna_mrna_promoters.gtf")

# counting promoter overlaps
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, type = "counts")

# row sum for each DBP
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)

# percent promoter overlap
num_peaks_df$percent_overlap <- round((num_peaks_df$peaks_overlapping_promoters/num_peaks_df$num_peaks*100),0)

write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/num_peaks_df.csv")

num_peaks_df
```

    ##     dbp num_peaks total_peak_length      ensembl_id      dbd  tf
    ## 1 CEBPZ      1892            740707 ENSG00000115816  Unknown Yes
    ## 2  CHD2     96965          49710914 ENSG00000173575 Myb/SANT  No
    ## 3  CTCF    316789         204671555 ENSG00000102974  C2H2 ZF Yes
    ## 4  ELF1    130614          83297027 ENSG00000120690      Ets Yes
    ## 5 EP300      6248           5503663 ENSG00000100393  Unknown  No
    ##   peaks_overlapping_promoters percent_overlap
    ## 1                        2343             124
    ## 2                       88968              92
    ## 3                       42328              13
    ## 4                      102344              78
    ## 5                         748              12

## results:
1) What can you determine from these overlaps? The number of the peaks
that overlap promoters for CEBPZ is greater than the number of peaks
itself. And the number of peaks that overlap promoters for CHD2 is
almost the same for the number that of peaks. I can say the DBPs likes
to bind near promoters.

## Now I want to compare the overlaps with lncRNA and mRNA promoters seperately

``` r
# last handy annotation will be lncRNA and mRNA gene IDs to subset
lncrna_gene_ids <- mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "lncRNA"]
mrna_gene_ids <-mrna_lncrna_genes$gene_id[mrna_lncrna_genes$gene_type == "protein_coding"]


#lncRNA promoter overlaps
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])

# mrna promoter overlaps
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])

num_peaks_df
```

    ##     dbp num_peaks total_peak_length      ensembl_id      dbd  tf
    ## 1 CEBPZ      1892            740707 ENSG00000115816  Unknown Yes
    ## 2  CHD2     96965          49710914 ENSG00000173575 Myb/SANT  No
    ## 3  CTCF    316789         204671555 ENSG00000102974  C2H2 ZF Yes
    ## 4  ELF1    130614          83297027 ENSG00000120690      Ets Yes
    ## 5 EP300      6248           5503663 ENSG00000100393  Unknown  No
    ##   peaks_overlapping_promoters percent_overlap
    ## 1                        2343             124
    ## 2                       88968              92
    ## 3                       42328              13
    ## 4                      102344              78
    ## 5                         748              12
    ##   peaks_overlapping_lncrna_promoters peaks_overlapping_mrna_promoters
    ## 1                                374                             1969
    ## 2                              18634                            70334
    ## 3                              14707                            27621
    ## 4                              21703                            80641
    ## 5                                429                              319

## results: 1) What is the difference in overlaps between mRNA and lncRNA promoters
Generally, the values of the peaks overlapping lncRNA promoters are much
smaller than that of mRNA promoters. Except for EP300 which has a higher
value of peaks overlapping lncRNA promoter to mRNA promoters.

## Now I am going to test if there is more binding over gene bodies than promoters
## I will seperate lncRNA and mRNA gene bodies to find the overlaps

``` r
# Note this takes 10-15 min
genebody_peak_counts <- count_peaks_per_feature(mrna_lncrna_genes, 
                                                consensus_peaks, 
                                                type = "counts")

# Now let's extract the overlaps the same way we did for promoters above

# All gene bodies
num_peaks_df$peaks_overlapping_genebody <- 
  rowSums(genebody_peak_counts)

# lncRNA gene bodies 
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])

# mRNA gene bodies
num_peaks_df$peaks_overlapping_mrna_genebody <- 
  rowSums(genebody_peak_counts[,mrna_gene_ids])

# let's take a look.
# cool lots of info let's save:
write_csv(num_peaks_df, "/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/num_peaks_df.csv")

num_peaks_df
```

    ##     dbp num_peaks total_peak_length      ensembl_id      dbd  tf
    ## 1 CEBPZ      1892            740707 ENSG00000115816  Unknown Yes
    ## 2  CHD2     96965          49710914 ENSG00000173575 Myb/SANT  No
    ## 3  CTCF    316789         204671555 ENSG00000102974  C2H2 ZF Yes
    ## 4  ELF1    130614          83297027 ENSG00000120690      Ets Yes
    ## 5 EP300      6248           5503663 ENSG00000100393  Unknown  No
    ##   peaks_overlapping_promoters percent_overlap
    ## 1                        2343             124
    ## 2                       88968              92
    ## 3                       42328              13
    ## 4                      102344              78
    ## 5                         748              12
    ##   peaks_overlapping_lncrna_promoters peaks_overlapping_mrna_promoters
    ## 1                                374                             1969
    ## 2                              18634                            70334
    ## 3                              14707                            27621
    ## 4                              21703                            80641
    ## 5                                429                              319
    ##   peaks_overlapping_genebody peaks_overlapping_lncrna_genebody
    ## 1                       2277                               341
    ## 2                     114983                             20647
    ## 3                     244288                             57948
    ## 4                     151987                             29953
    ## 5                       4829                              1562
    ##   peaks_overlapping_mrna_genebody
    ## 1                            1936
    ## 2                           94336
    ## 3                          186340
    ## 4                          122034
    ## 5                            3267

## results: 1) Do my proteins have more overlaps with promoters or genebodies?
Generally, there are more peaks overlapping genebodies than promoters.
When lncRNA and mRNA are separated, the values of the peaks overlapping
the promoters are smaller than the values of the peaks overlapping the
genebody except for CEBPZ.

## It is nice and all to find overlaps, but I am interested in how many proteins
## bind a specific promoter. I will use my handy “occurence” parameter in
## " count peaks per feature"

``` r
#using the type=occurrence instead of type=count
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, consensus_peaks, 
                                               type = "occurrence")

# Let's double check that all lncrna & mrna genes are accounted for:
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

# Great we will use this quite a bit moving forward so let's write it out! 
write.table(promoter_peak_occurence, "my_results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# First make sure promoter_peak_occurrence and lncrna_mrna_promoters are in the same order
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))

#how many proteins are bound by a promoter.
# Now let's use the 'data.frame()' fucntion. Set up a bunch of colnames and populate them.
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))


# Let's write out this data frame -- all this code & run time is now encapsulated in one .csv :)
write_csv(peak_occurence_df, "my_results/peak_occurence_dataframe.csv")
table(peak_occurence_df$number_of_dbp)
```

    ## 
    ##     0     1     2     3     4 
    ## 24271  5895  5363  1270    15

## results: I find the max number of proteins on a promoter to be 4


## Now I want to start plotting my results
## First I will see if there is a realtionship between peak number and total DNA covered

``` r
# let's make this simple plot first: number of peaks -vs- total peak length
num_peaks_df <- read_csv('/scratch/Shares/rinnclass/CLASS_2023/Christopher/CLASS_2023/CLASSES/num_peaks_df.csv')
```

    ## Rows: 5 Columns: 13
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (4): dbp, ensembl_id, dbd, tf
    ## dbl (9): num_peaks, total_peak_length, peaks_overlapping_promoters, percent_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
ggplot(num_peaks_df, aes(x = num_peaks, 
                         y = total_peak_length,
                         color = dbp)) +
  geom_point() +
  geom_text(aes(label = dbp, color = dbp), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5) +
  scale_color_discrete()
```

![](christopher_new_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# We have a result! The more peaks the more genome space is covered.
```

## Now I want to color my plot by whether the protein is a TF or not.

``` r
# Let's plot the num_peaks -vs- coverage for Zinc finger proteins (C2H2 ZF)
ggplot(num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                 color = dbd == "C2H2 ZF")) +
  geom_point() +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5) 
```

![](christopher_new_files/figure-gfm/none%20of%20my%20proteins%20are%20TF-1.png)<!-- -->

``` r
# Let's plot the num_peaks -vs- coverage for Ets
ggplot(num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                 color = dbd == "Ets")) +
  geom_point() +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5) 
```

![](christopher_new_files/figure-gfm/none%20of%20my%20proteins%20are%20TF-2.png)<!-- -->

``` r
# Let's plot the num_peaks -vs- coverage for Unknown
ggplot(num_peaks_df, aes(x = num_peaks, 
                 y = total_peak_length,
                 color = dbd == "Unknown")) +
  geom_point()  +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5) 
```

![](christopher_new_files/figure-gfm/none%20of%20my%20proteins%20are%20TF-3.png)<!-- -->

## I want to make a histogram of the number of peaks for each of my proteins

``` r
# frequencty of peak numbers across DBPs

ggplot(num_peaks_df, aes(x = num_peaks)) +
  geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](christopher_new_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# Result: peak number is represented by right gaussian.
# let's add more bins

ggplot(num_peaks_df, aes(x = num_peaks, fill = tf)) +
   geom_histogram(bins = 5, position = "dodge")
```

![](christopher_new_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

## Now I want to facet this by the type of DNA binding domain my protein has.

``` r
# Now let's set dbds to these three DBDs.
dbds <- c("C2H2 ZF", "Ets", "Unknown", "Myb/SANT") 

# Now the plotting beings and good example of running functions in aes()
ggplot(num_peaks_df %>% 
         filter(dbd %in% dbds),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_grid(dbd ~ .) + 
  geom_point()  +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.2), 
            vjust = 1.4, hjust = 1) 
```

![](christopher_new_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#We can also facet grid another parameter -- let's say tf
ggplot(num_peaks_df %>% 
         filter(dbd %in% dbds),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_grid(dbd ~ tf) + 
  geom_point() +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.2), 
            vjust = 0.6, hjust = 0.5)
```

![](christopher_new_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
# Here we said plot DBD as rows and tf status as columns -- 
# We can easily flip this too:


ggplot(num_peaks_df %>% filter(dbd %in% dbds),
       aes(x = num_peaks, y = total_peak_length )) +
  facet_grid(tf ~ dbd) + 
  geom_point() +
  geom_text(aes(label = dbp), 
            position = position_dodge(width = 0.2), 
            vjust = 0.6, hjust = 0.5)
```

![](christopher_new_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->
