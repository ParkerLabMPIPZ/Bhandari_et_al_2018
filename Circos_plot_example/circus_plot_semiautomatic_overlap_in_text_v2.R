## circus plot for looking at the overlap between given sets of genes ###

## adjusted on 13.07.2018
## author Dmitry + help pages from circlize
## Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014
## https://jokergoo.github.io/circlize_book/book/

rm(list = ls())
## libraries
library(circlize)
library(colorspace)


###### Known issues #####
# I have noticed that circlize fails to work with very large datasets of >10000 genes as well as with very small datasets (~20 genes).

# you will notice a warning message ~ "Note: 20 points are out of plotting region in sector ... "
# this is "normal" behaviour for circlize package. Only a couple of points are not drawn on the plot but
# this is not of a concern for large sets

# colorspace package is used to produce color palettes and not well-recognised Rcolorbrewer because
# the later is suitable for small number of separate colours (up to 12)

############################################
###### REQUIRED INPUT AND ADJUSTMENT #######
############################################
# IMPORTANT
# files with gene sets should contain at least one column with the name "geneID", which contains Arabidopsis AGI
# codes in the format AT1G11111
# input format - tab delimited


# put all geneset files (one file per set) into a folder "genesets_for_overlap" inside your working directory
folder_with_genesets <- "./genesets_for_overlap"

# provide exact name of the file containing names of genes in the main set of interest 
name_of_file_with_main_dataset <- "Cui_EP_estradiol_12_24h.txt"

# provide name of the output graphics in pdf format
output_circos_plot_file_name <- "test.pdf"

# create output folder to write files with overlaps between datasets
# taken from https://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
mainDir <- getwd()
subDir <- "output_overlap_between_sets"
dir.create(file.path(mainDir, subDir))

#############################################################
#### THE REST DOES NOT NORMALLY REQUIRE MODIFICATION ########
#############################################################

#### read user-provided sets of genes for testing overlap between them
geneset_file_names <- list.files(folder_with_genesets)
geneset_names <- unlist(strsplit(geneset_file_names,"[.]"))[seq(from = 1, to = length(geneset_file_names)*2, by=2)]
genesets <- list()

for (i in 1:length(geneset_file_names)){
  path_to_file <- paste(folder_with_genesets,"/",geneset_file_names[i],sep = "")
  genesets[i] <- read.delim(path_to_file,header=TRUE)
}


#### prepare dataframe for circos plot
df_for_circos<-data.frame()

for (i in 1:length(geneset_names)){
  genes_in_given_dataset <- unlist(genesets[i])
  nr_genes_in_given_dataset <- length(genes_in_given_dataset)
  temp <- data.frame(gene_set_name = rep(geneset_names[i],nr_genes_in_given_dataset),
                     GeneID = genes_in_given_dataset,
                     integer_representation = 1:nr_genes_in_given_dataset)
  df_for_circos <- rbind(df_for_circos,temp)
  rm(temp)
}

sets_for_overlap <- levels(df_for_circos$gene_set_name)

#### set colours for sectors on the circos plot
length(sets_for_overlap) #brewer pallete has limited number of colours -> use colorspace package

colour_for_sectors <- sequential_hcl(length(sets_for_overlap))

#### set colours for sectors on the circos plot
colour_for_links <- diverge_hcl(length(sets_for_overlap))

#### circos plot itself
df_for_circos$y <- 0.5 ## y axis is required to initialize layout, just to give a position


circos.clear()

pdf(output_circos_plot_file_name, width = 20, height = 20)

circos.par("track.height" = 0.08) # set global parameters with circos.par()
circos.initialize(factors = df_for_circos$gene_set_name, x = df_for_circos$integer_representation) # initialize a circular layout - factor for sectors

#initialize track
circos.track(factors = df_for_circos$gene_set_name, ylim = c(0.48, 0.52),
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(4, "mm"), 
                           CELL_META$sector.index,
                           cex=1.5,
                           facing="bending.inside")
               circos.axis(h = "bottom",labels.cex = 0.6)
             })



circos.trackPoints(df_for_circos$gene_set_name, df_for_circos$integer_representation, df_for_circos$y, col = colour_for_sectors, pch = 16, cex = 2)


# add links for the overlap
name_main_geneset_to_overlap <- unlist(strsplit(name_of_file_with_main_dataset,"[.]"))[1]
main_geneset_to_overlap <- df_for_circos$GeneID[df_for_circos$gene_set_name==name_main_geneset_to_overlap]
sets_for_overlap_without_main <- which(sets_for_overlap!=name_main_geneset_to_overlap)


for (k in sets_for_overlap_without_main){
  intersect_clusters<-intersect(df_for_circos$GeneID[df_for_circos$gene_set_name==sets_for_overlap[k]],main_geneset_to_overlap)
  if (length(intersect_clusters)>0) {
    for (i in 1:length(intersect_clusters)){
      circos.link(sets_for_overlap[k], df_for_circos$integer_representation[df_for_circos$gene_set_name==sets_for_overlap[k] & 
                                                                   df_for_circos$GeneID==intersect_clusters[i]],
                  name_main_geneset_to_overlap, df_for_circos$integer_representation[df_for_circos$gene_set_name==name_main_geneset_to_overlap & 
                                                                         df_for_circos$GeneID==intersect_clusters[i]],
                  col=colour_for_links[k],
                  lwd = 1)
    }
  }
  write.table(intersect_clusters,
              file.path(mainDir, subDir, paste("genes_in_",name_main_geneset_to_overlap,"_and_",sets_for_overlap[k],".txt",sep="")),
              col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  rm(intersect_clusters)
}

dev.off()
