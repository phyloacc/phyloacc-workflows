library(tidyverse)
library(grid)

#############################################################################

chr_list = str_c("chr", as.character(1:22))

chr_info_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/00-human-ref-ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
conserved_sites_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/04-phylop-bedgraph/241-mammalian-2020v2b.phylop-scores.fdr-adj-0.05.conserved.tsv"

chr_info = read_tsv(chr_info_file, col_names=c("chr", "length", "location", "line.length", "line.length.nl")) %>% 
  select(chr, length) %>%
  mutate(chr = paste0("chr", chr)) %>%
  filter(chr %in% chr_list)

conserved_sites = read_tsv(conserved_sites_file, col_names=c("chr", "start", "end", "sig.code")) %>% select(chr, start)

for(chrom in chr_list){
  conserved_sites_chr = conserved_sites %>% 
    filter(chr==chrom)
  
  null_model = runif(nrow(conserved_sites_chr), min=1, max=chr_info[chr_info$chr==chrom,]$length)
  
  # Perform the Kolmogorov-Smirnov test
  ks_test = ks.test(conserved_sites_chr$start, null_model)
  
  corrected_pval = ks_test$p.value / length(chr_list)

  outline = paste(chrom, ks_test$p.value, corrected_pval)
  
  if(corrected_pval < 0.01){
    outline = paste(outline, "*")
  }
  
  print(outline)
  
  conserved_elements_file_all = paste0("/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/07-conserved-elements-chr/", chrom, "/80/241-mammalian-2020v2b.phylop-conserved-windows.", chrom, ".0.05.80.40.bed")
  conserved_elements_file_241 = paste0("/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/08-conserved-elements-chr-fasta/", chrom, "/241-mammalian-2020v2b.phylop-conserved-windows.", chrom, ".0.05.80.40.241.bed")
  conserved_elements = read_tsv(conserved_elements_file_241, col_names=c("chr", "start", "end"))
  
  chr_len = chr_info[chr_info$chr==chrom,]$length
  site_win_size = 100000
  elem_win_size = 10000000
  x_interval = 10000000
  p = ggplot(conserved_sites_chr, aes(x=start)) +
    geom_histogram(binwidth=site_win_size, fill="#db6d00") +
    geom_histogram(data = conserved_elements, aes(x = start, y = -..count..), binwidth = elem_win_size, fill = "#009292", color = "#000000", size=0.25) +
    scale_x_continuous(limits = c(0, chr_len), breaks = seq(0, chr_len, by = x_interval),
                       labels = function(x) paste0(x/1000000)) +
    scale_y_continuous(breaks = c(seq(-1500, 0, by = 1000), seq(0, 21000, by = 5000)),
                       labels = function(x) abs(x)) +
    coord_cartesian(ylim = c(-1500, 21000)) +
    geom_hline(yintercept=0) +
    xlab("Position (Mb)") +
    ylab("Count") +
    #labs(y = "# Conserved sites per 100kb window") +
    annotate("text", x = 60000000, y = 20750, label = "Conserved sites per 100kb window") +
    annotate("text", x = 65000000, y = -1500, label = "Conserved elements per 10Mb window") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  print(p)
  
  outfile = paste0("../figs/cons-dist-", chrom, ".pdf")
  ggsave(filename=outfile, plot=p, height=4, width=6, units="in")
  stop("OK")
}




# grid.text("# Conserved elements\nper 10Mb window", x = unit(0, "npc") + unit(0.5, "lines"), 
#           y = unit(0.2, "npc"), rot = 90)

