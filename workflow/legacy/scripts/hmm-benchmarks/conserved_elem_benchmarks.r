
library(tidyverse)
#library(plotly)
source("/n/home07/gthomas/env/pylib/core/r/design.r")

hmm_data_file = "combined.log"
hmm_data = read_delim(hmm_data_file, delim="\t")

hmm_data = hmm_data %>% mutate(t.rat = t0_0 / t1_1) %>%
  mutate(t.diff = t0_0 - t1_1) %>%
  mutate(t.prod = t0_0 * t1_1) %>%
  mutate(t.str = paste(t0_0, t1_1, sep="-")) %>%
  mutate(e.rat = e0_0 / e1_1) %>%
  mutate(e.diff = e0_0 - e1_1) %>%
  mutate(e.prod = e0_0 * e1_1) %>%
  mutate(e.str = paste(e0_0, e1_1, sep="-"))


# Get a list of all .bed files in the working directory
# hmm_bed_files <- list.files(pattern = "\\.bed$")


# Loop over the file list
# for (bed_file in hmm_bed_files) {
#   print(bed_file)
#   if (file.size(bed_file) == 0) next
#   
#   parts = unlist(strsplit(bed_file, "-"))
#   cur_t_str = paste(parts[4], parts[5], sep = "-")
#   cur_e_str = paste(parts[6], parts[7], sep = "-")
#   
#   cur_bed_data = read_delim(bed_file, delim="\t", col_names=c("chr", "start", "end", "length", "cons.sites"))
#   cur_bed_data = cur_bed_data %>% mutate(perc.cons.sites = cons.sites / length)
#   
#   cons_sites_in_cons_elem = sum(cur_bed_data$cons.sites)
#   avg_perc_cons_sites_per_elem = mean(cur_bed_data$perc.cons.sites)
#   
#   hmm_data = hmm_data %>% 
#     mutate(cons.sites.in.cons.elem = ifelse(t.str == cur_t_str & e.str == cur_e_str, cons_sites_in_cons_elem, cons.sites.in.cons.elem)) %>%
#     mutate(avg.perc.cons.sites.per.elem = ifelse(t.str == cur_t_str & e.str == cur_e_str, avg_perc_cons_sites_per_elem, avg.perc.cons.sites.per.elem))
# }
# 
hmm_data = hmm_data %>% mutate(perc.cons.sites.in.cons.elem = cons.sites.in.cons.elems / conserved.sites)

hmm_data_t = hmm_data %>% filter(e0_0 == 0.8 & e1_1 == 0.8)
hmm_data_e = hmm_data %>% filter(t0_0 == 0.5 & t1_1 == 0.5)


p = ggplot(hmm_data_t, aes(perc.cons.sites.in.cons.elem, avg.perc.cons.sites.per.elem, color=avg.cons.elem.len)) +
  geom_point() +
  bartheme()
print(p)

p = ggplot(hmm_data_t, aes(perc.cons.sites.in.cons.elem, avg.perc.cons.sites.per.elem, color=avg.cons.elem.len)) +
  geom_point() +
  bartheme()
print(p)
             
stop("OK")


hmm_data %>% filter(avg.inter.elem.dist < 5e5) %>%
  ggplot(aes(t.str, e.str, fill=avg.inter.elem.dist)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    xlab("Same state transition probs") +
    ylab("Same state emission probs") +
    labs(fill="Avg. inter-element\ndistance") +
    bartheme() +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1, size=4),
          axis.text.y = element_text(size=4),
          legend.title=element_text())


hmm_data %>% filter(avg.cons.elem.len < 1000 & avg.inter.elem.dist < 1000) %>%
  ggplot(aes(avg.cons.elem.len, avg.inter.elem.dist, color=final.cons.elems)) +
    geom_point() +
    xlab("Avg. conserved element length") +
    ylab("Avg. inter-element distance") +
    labs(color="# conserved elements") +
    bartheme() +
    theme(legend.title=element_text())


hmm_data_t %>% filter(final.cons.elems < 900000) %>%
  ggplot(aes(t0_0, t1_1, color=final.cons.elems)) +
  geom_point() +
  bartheme()


#plot_ly(hmm_data, x = ~t0_0, y = ~t1_1, z = ~e0_0, color = ~e1_1, size = ~final.cons.elems, type = "scatter3d", mode = "markers")


# plot transition probs with % of cons sites in cons elems and % of each cons elem that is cons
# clustering in genomes/genes
# clustering in a series
