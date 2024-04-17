library(tidyverse)

#############################################################################

chr_list = str_c("chr", as.character(1:22))
#chr_list = str_c("chr", as.character(20:22))

conserved_sites_file = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/04-phylop-bedgraph/241-mammalian-2020v2b.phylop-scores.fdr-adj-0.05.conserved.tsv"
depth_dir = "/n/holylfs05/LABS/informatics/Users/gthomas/phyloacc-workflows/data/01-zoonomia-aln/03-aln-depth/"
all_depth_files = list.files(depth_dir)

conserved_sites = read_tsv(conserved_sites_file, col_names=c("chr", "start", "end", "sig.code")) %>% select(chr, start)

plot_list = list()

for(chrom in chr_list){
  print(chrom)
  cur_depth_files = str_c(depth_dir, all_depth_files[str_detect(all_depth_files, paste0("^", chrom, "\\."))])

  print("# reading files...")
  cur_depth_sites = cur_depth_files %>%
    map_df(~{
      print(.x)  # This will print the name of the file
      read_tsv(.x, col_names=c("chr", "start", "end", "depth")) %>%
        select(chr, start, depth)
    })
  
  print("# subsetting conserved sites...")
  conserved_sites_chr = conserved_sites %>% 
    filter(chr==chrom) %>%
    mutate(conserved = TRUE)
  
  print("# combining labels...")
  cur_depth_sites = cur_depth_sites %>% 
    left_join(conserved_sites_chr, by = c("chr", "start")) %>% 
    mutate(conserved = replace_na(conserved, FALSE))
  
  print("# plotting...")
  p = ggplot(cur_depth_sites, aes(x=conserved, y=depth, group=conserved)) +
    geom_boxplot(outlier.shape=NA, fill="#cccccc") +
    xlab("Conserved?") +
    ylab("Alignment depth") +
    ggtitle(chrom) +
    theme_classic()
  #print(p)
  
  outfile = paste0("../figs/cons-v-depth-", chrom, ".pdf")
  ggsave(filename=outfile, plot=p, height=4, width=4, units="in")
  
  stop("# OK!")
}

#############################################################################