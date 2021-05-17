#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
})



check_and_collapse = function(raw_DNA_folder, raw_RNA_folder, output.name="qc_output") {
  
  ### glossary:
  # read - numbers in column 7, usually number of raw sequencing reads
  # unit - column 8 (if provided, currently ignored), UMIs/regions/barcodes
  # element - column 4 (aka sequence, oligo, insert, tile...)
  
  # find files
  if(missing("output.name")) stop("Argument 'output.name' is missing!")
  raw_DNA_files = list.files(raw_DNA_folder, pattern = '.bed(.gz)?',full.names = T)
  raw_RNA_files = list.files(raw_RNA_folder, pattern = '.bed(.gz)?',full.names = T)
  
  ### check number of samples
  DNA_len = length(raw_DNA_files)
  RNA_len = length(raw_RNA_files)
  if(RNA_len==0 | DNA_len ==0) stop("Couldn't find DNA or RNA files (BED or BED.GZ)")
  if(RNA_len==DNA_len) {
    message("Found ", DNA_len, " DNA and RNA replicates.") } else {
      if(RNA_len>DNA_len & DNA_len==1) {
        message("Found only one DNA reference for ", RNA_len, " RNA replicates.") } else {
          stop("Found ", list.files(raw_DNA_folder, pattern = '.bed(.gz)?'),
               " DNA samples and ", list.files(raw_RNA_folder, pattern = '.bed(.gz)?'), 
               " RNA samples. That's unacceptable, aborting.")
        }
    }
  
  ### make output folder
  output.dir = file.path("report",output.name)
  dir.create(output.dir,showWarnings = F,recursive = T)
  dir.create(paste0(output.dir,"/","DNA"),showWarnings = F,recursive = T)
  dir.create(paste0(output.dir,"/","RNA"),showWarnings = F,recursive = T)
  # round(file.size(raw_DNA_files[replicate])/1024^2,1) # file size in MB, if needed
  
  ### collapse (if need be)
  raw_files = c(raw_DNA_files,raw_RNA_files)
  collapsed_DNA = paste0(output.dir,sub(".*(\\/.*\\/.*)$","\\1",
                                        gsub("\\/\\/","\\/",raw_DNA_files)))
  collapsed_RNA = paste0(output.dir,sub(".*(\\/.*\\/.*)$","\\1",
                                        gsub("\\/\\/","\\/",raw_RNA_files)))
  collapsed_files = c(collapsed_DNA,collapsed_RNA)
  
  
  
  for(i in seq_along(raw_files)) {
    
    if(file.exists(collapsed_files[i])) next()
    message("Processing: ", raw_files[i])
    header = grepl("start|end|strand", readLines(raw_files[i],n = 1))
    
    system(paste0(ifelse(grepl("\\.gz$", raw_files[i]), "gzcat ", "cat "), 
                  raw_files[i], 
                  ifelse(header," | tail -n+2", ""),
                  " | sort -k4 ", 
                  " | bedtools groupby -full -g 4 -c 7 -o sum,count",
                  " | cut -f1-4,9,10 ",
                  " | gzip > ",
                  collapsed_files[i]))
    
  }
  
  out = list(collapsed_DNA,collapsed_RNA)
  names(out) = c("DNA","RNA")
  return(out)
}


mpra_wrap = function(raw_DNA_folder, raw_RNA_folder, threshold = -Inf, output.name) {
  
  # check file numbers, collapse to element (col 4) level, 
  # output paths to collapsed files
  collapsed_files = check_and_collapse(raw_DNA_folder = raw_DNA_folder,
                                       raw_RNA_folder = raw_RNA_folder,
                                       output.name = output.name)
  output.dir = file.path("report",output.name)
  
  ##########################
  ###### DNA section #######
  ##########################
  

  threshold=as.numeric(threshold)
  recalculate_threshold = is.infinite(threshold)
  dna_filter = list(); threshold.reads=list(); threshold.molecules=list()
  # message(recalculate_threshold, " recalculate thr")
  
  for(replicate in seq_along(collapsed_files$DNA)) {
    
    dna = read_tsv(collapsed_files$DNA[replicate],
                   col_names = c('chr','start','end','name',
                                 'DNA.reads','DNA.molecules'), 
                   col_types = "ciicii")
    already_collapsed_dna = !any(dna$DNA.molecules>1)
    
    ### thresholding: currently only total DNA reads, not units (ie UMI, regions etc)
    # if no threshold provided, use Q1-1.5*IQR
    if(recalculate_threshold) {
      threshold.reads[[replicate]] = as.integer(10^(quantile(log10(dna$DNA.reads), probs = c(0.25)) - 
                                                      1.5*IQR(log10(dna$DNA.reads)))) 
    } else {
      threshold.reads[[replicate]] = threshold
    }
    if(!already_collapsed_dna) {
      threshold.molecules[[replicate]] = as.integer(10^(quantile(log10(dna$DNA.molecules), probs = c(0.25)) - 
                                                          1.5*IQR(log10(dna$DNA.molecules))))
    }
    
    # message(replicate)
    # message(threshold.reads[[replicate]])
    
    ### histograms of DNA
    
    hist.reads = dna %>%
      ggplot(aes(DNA.reads))+
      geom_histogram(bins=50) +
      geom_vline(xintercept = threshold.reads[[replicate]], color="red") +
      scale_x_log10() +
      labs(title = paste0("DNA replicate ", replicate))+
      xlab('DNA reads per element')
    
    if(!already_collapsed_dna) {
      
      hist.molecules = dna %>%
        ggplot(aes(DNA.molecules))+
        geom_histogram(bins=50) +
        geom_vline(xintercept = threshold.molecules[[replicate]], color="red") +
        # facet_grid(~key) +
        scale_x_log10() +
        # labs(title = paste0("DNA replicate ", replicate))+
        xlab('DNA molecules per element')
      
      hist.read.per.molecule = dna %>%
        mutate(read.per.molecule = DNA.reads/DNA.molecules) %>% 
        ggplot(aes(read.per.molecule))+
        geom_histogram(bins=50) +
        scale_x_log10() +
        xlab('DNA reads per molecule')
      
      # save plot & filter
      plot.hist.dna = plot_grid(hist.reads, hist.molecules, hist.read.per.molecule, nrow = 1,align = 'h')
      ggsave(plot = plot.hist.dna, 
             file.path(output.dir, paste0("DNA_rep",replicate,"_hist_plot.png")),
             w=8,h=2)
      
      dna_filter[[replicate]] = dna %>% 
        mutate(keep = DNA.reads >= threshold.reads[[replicate]] &
                 DNA.molecules >= threshold.molecules[[replicate]]) 
      
    } else {
      
      # save plot & filter
      ggsave(plot = hist.reads, 
             file.path(output.dir, paste0("DNA_rep",replicate,"_hist_plot.png")),
             w=4,h=2)
      
      dna_filter[[replicate]] = dna %>% 
        mutate(keep = DNA.reads >= threshold.reads[[replicate]]) 
    }
    
  }
  
  ##########################
  ###### RNA section #######
  ##########################
  
  rna_filter = list()
  
  for(replicate in seq_along(collapsed_files$RNA)) {
    
    rna = read_tsv(collapsed_files$RNA[replicate],
                   col_names = c('chr','start','end','name',
                                 'RNA.reads','RNA.molecules'), 
                   col_types = "ciicii")
    
    already_collapsed_rna = !any(rna$RNA.molecules>1)
    
    ### histograms of RNA
    
    hist.reads = rna %>%
      ggplot(aes(RNA.reads))+
      geom_histogram(bins=50) +
      scale_x_log10() +
      labs(title = paste0("RNA replicate ", replicate))+
      xlab('RNA reads per element')
    
    if(!already_collapsed_rna) {
      
      hist.molecules = rna %>%
        ggplot(aes(RNA.molecules))+
        geom_histogram(bins=50) +
        scale_x_log10() +
        xlab('RNA molecules per element')
      
      hist.read.per.molecule = rna %>%
        mutate(read.per.molecule = RNA.reads/RNA.molecules) %>% 
        ggplot(aes(read.per.molecule))+
        geom_histogram(bins=50) +
        scale_x_log10() +
        xlab('RNA reads per molecule')
      
      # save plot & filter
      plot.hist.rna = plot_grid(hist.reads, hist.molecules, hist.read.per.molecule, nrow = 1,align = 'h')
      ggsave(plot = plot.hist.rna, 
             file.path(output.dir, paste0("RNA_rep",replicate,"_hist_plot.png")),
             w=8,h=2)
      
      rna_filter[[replicate]] = rna
      
    } else {
      
      # save plot & filter
      ggsave(plot = hist.reads, 
             file.path(output.dir, paste0("RNA_rep",replicate,"_hist_plot.png")),
             w=4,h=2)
      
      rna_filter[[replicate]] = rna
    }
  }
  
  
  
  ##########################
  ##### ratios section #####
  ##########################
  
  
  
  stats = list(); ratios=list()
  
  # situation where the number of RNA replicates != DNA replicates
  # is allowed, but not implemented here
  for(replicate in seq_along(collapsed_files$RNA)) {
    
    message(replicate)
    current_rna = rna_filter[[replicate]]
    current_dna = dna_filter[[replicate]]
    current_dna_filter = current_dna %>% 
      filter(keep) %>% 
      mutate(reads.per.molecule = DNA.reads/DNA.molecules)
    
    ### calculate ratios (assuming "DNA is right", ie left_join)
    ratios[[replicate]] = current_dna_filter %>%
      left_join(current_rna,by="name") %>%
      mutate(DNA.norm = log2( (DNA.reads+1)/sum(DNA.reads)*1e6 ),
             RNA.norm = log2( (RNA.reads+1)/sum(RNA.reads)*1e6 ),
             log.ratio = RNA.norm - DNA.norm,
             replicate)
    
    
    
    stats[[replicate]] = tibble(
      
      output.name, 
      replicate, 
      
      dna_read_threshold = threshold.reads[[replicate]],
      dna_molecule_threshold = ifelse(!already_collapsed_dna, threshold.reads[[replicate]], NA),
      
      saturation_dna = ifelse(!already_collapsed_dna, 
                              round(1-sum(current_dna$DNA.molecules)/sum(current_dna$DNA.reads),2),NA), 
      
      dna_elements = nrow(current_dna),
      dna_elements_filter = sum(current_dna$keep),
      
      dna_molecules = ifelse(!already_collapsed_dna, 
                             sum(current_dna$DNA.molecules),NA),
      dna_molecules_filter = ifelse(!already_collapsed_dna, 
                                    sum(current_dna_filter$DNA.molecules),NA),
      
      dna_reads = sum(current_dna$DNA.reads),
      dna_reads_filter = sum(current_dna_filter$DNA.reads),
      
      median_dna_molecules_per_elem_filter = ifelse(!already_collapsed_dna, 
                                                    median(current_dna_filter$DNA.molecules),NA),
      median_dna_reads_per_molecule_filter = ifelse(!already_collapsed_dna, 
                                                    round(median(current_dna_filter$reads.per.molecule),2),NA),
      
      saturation_rna = ifelse(!already_collapsed_rna, 
                              round(1-sum(current_rna$RNA.molecules)/sum(current_rna$RNA.reads),2),NA), 
      rna_reads = sum(current_rna$RNA.reads),
      rna_molecules = ifelse(!already_collapsed_rna, 
                             sum(current_rna$RNA.molecules),NA)
    )
    
  }
  
  
  # write result to files
  write_tsv(bind_rows(stats),file.path(output.dir,"complete_stats.tsv"))
  complete_ratios = bind_rows(ratios) %>% 
    mutate(label = paste0("rep",replicate))
  write_tsv(complete_ratios,gzfile(file.path(output.dir,"complete_ratios.tsv.gz")))
  
  
  ### make correlation plots (RNA,DNA,ratios)
  # suppressing warnings because of pairwise incomplete observations
  pl.dna = complete_ratios %>% 
    select(label,DNA.norm,name) %>% 
    spread(label,DNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "DNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_DNA.png"),pl.dna,w=7,h=7))
  
  pl.rna = complete_ratios %>% 
    select(label,RNA.norm,name) %>% 
    spread(label,RNA.norm) %>% 
    select(-name) %>% 
    ggpairs(title = "RNA RPM correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_RNA.png"),pl.rna,w=7,h=7))
  
  pl.ratios = complete_ratios %>% 
    select(label,log.ratio,name) %>% 
    spread(label,log.ratio) %>% 
    select(-name) %>% 
    ggpairs(title = "ratio correlation")
  suppressWarnings(ggsave(file.path(output.dir, "cor_ratio.png"),pl.ratios,w=7,h=7))
  
}  


# Parse command-line arguments
parser <- ArgumentParser()
parser$add_argument("-d", "--dnas", required=T, help="Folder with Raw DNA counts")
parser$add_argument("-r", "--rnas", required=T, help="Folder with Raw RNA counts")
# parser$add_argument("--method", required=F, choices=c("MPRA", "STARR"),  
#                     default="MPRA", help="Specify the experiment type.")
parser$add_argument("-th","--threshold-reads", required=F,   
                    default="-Inf", help="Threshold used to filter out low DNA count elements.")
parser$add_argument("-o", "--outname", required=T, 
                    help="Output rootname for newly generated QC files.")

args <- parser$parse_args()

# load libraries after successfully parsing the command-line arguments
suppressPackageStartupMessages({
  library(tidyverse)
  library(GGally)
  library(cowplot)
})
theme_set(theme_bw())
mpra_wrap(raw_DNA_folder = args$dnas, 
          raw_RNA_folder = args$rnas, 
          threshold=args$threshold_reads, 
          # method=args$method, 
          output.name = args$outname)