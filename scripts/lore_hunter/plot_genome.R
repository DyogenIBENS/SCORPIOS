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

if (opt$ncolors == 1){
	colors <- c('#FF0000')
}

if (opt$ncolors == 3){
	colors <-  c('#FF0000', "#332288", "#44AA99") #"#DDCC77", "#CC6677", "#117733", "#88CCEE", "#AA4499")
}

if (opt$ncolors == 2){
	colors <-  c('#FF0000', "#332288") #"#DDCC77", "#CC6677", "#117733", "#88CCEE", "#AA4499")
}


ideogram(karyotype = karyotype, overlaid = overlay_data, output = opt$outfile, width = 180, colorset1=colors)
         # colorset1= colors)
