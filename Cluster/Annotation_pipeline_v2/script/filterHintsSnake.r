#!/usr/local/bioinfo/R/3.2.2/bin/Rscript --vanilla
# -*- coding: utf-8 -*-

library("optparse")
library("stringr")


## define options
option_list = list(
  make_option(c("-s", "--specie"), type="character", default=NULL, help="input species name", metavar="filename"),
  make_option(c("-p", "--path"), type="character", default=NULL, help="path of working dir", metavar="path"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="path of output directory", metavar="output")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$specie)){
  print_help(opt_parser)
  stop("At least one argument must be supplied specie.\n", call.=FALSE)
}
if (is.null(opt$path)){
  print_help(opt_parser)
  stop("At least one argument must be supplied path.\n", call.=FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("At least one argument must be supplied output.\n", call.=FALSE)
}

path = opt$path
name = opt$specie
output = opt$out

setwd(path)


# 70-15
hints.raw <- read.table(pipe(paste("awk '/mult=/' hints_",name,".raw.bam", sep = "")), stringsAsFactors=F)
hints.eval <- as.data.frame(str_split_fixed(hints.raw$V9, ";", 3))

# SET minimum of reads required for an intron hint
min.reads <- 20

hints.eval.mult <- as.numeric(gsub("mult=","",hints.eval$V1))
hints.filtered <- hints.raw[hints.eval.mult >= min.reads,]

write.table(hints.filtered, output, sep="\t", quote=F, row.names=F, col.names=F)
