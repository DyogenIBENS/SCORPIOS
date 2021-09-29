library(RIdeogram)
library(optparse)

option_list = list(
  make_option(c("-k", "--karyo.file"), type="character", help="karyotype file for RIdeogram"),
  make_option(c("-f", "--features"), type="character", help="RIdeogram overlaid data file."),
  make_option(c("-o", "--outfile"), type="character", default="out_karyotype.svg",
              help="Output figure filename. [default= %default]"),
  make_option(c("-c", "--ncolors"), type="integer", default=1,
              help="Number of colors to use (= number of features to plot.) [default= %default]")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

karyotype <- read.table(opt$karyo.file, sep = "\t", header = T, stringsAsFactors = F)
overlay_data <- read.table(opt$features, sep = "\t", header = T, stringsAsFactors = F)

palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f",
             "#bcbd22","#17becf")

if (opt$ncolors == 1){
	colors <- c('#FF0000')
} else if (opt$ncolors > 1 && opt$ncolors < 11 ){
  colors <- palette[1:opt$ncolors]
} else {
  write("Error: can only color karyotype with up to 10 different colors", stderr())
  quit(status=11)
}

# if (opt$ncolors == 3){
# 	colors <-  c('#FF0000', "#332288", "#44AA99") #"#DDCC77", "#CC6677", "#117733", "#88CCEE", "#AA4499")
# }

# if (opt$ncolors == 2){
# 	colors <-  c('#FF0000', "#332288") #"#DDCC77", "#CC6677", "#117733", "#88CCEE", "#AA4499")
# }


ideogram(karyotype = karyotype, overlaid = overlay_data, output = opt$outfile, colorset1=colors)
