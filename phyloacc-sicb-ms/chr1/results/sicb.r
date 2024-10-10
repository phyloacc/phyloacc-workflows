
library(tidyverse)
here::i_am("results/sicb.r")
library(here)
library(viridis)
source(here("results", "lib", "design.r"))

bf1_thresh = 5
bf2_thresh = 5

max_bf1_disp = 200
min_bf2_disp = -500

results_raw = read_tsv(here("results", "elem_lik.txt"), comment="#")
results = results_raw %>% 
  mutate(above.bf1 = logbf1 > bf1_thresh) %>%
  mutate(above.bf2 = logbf2 > bf2_thresh) %>%
  mutate(accelerated = ifelse(logbf1 > bf1_thresh & logbf2 > 5, "Accelerated", "Not accelerated"))
  


p = results %>%
  filter(logbf1 < max_bf1_disp & logbf2 > min_bf2_disp) %>%
  ggplot(aes(x=logbf1, y=logbf2, color=accelerated)) +
  geom_point(alpha=0.3) +
  geom_vline(xintercept=bf1_thresh, color="#666666", linetype="dashed") + 
  geom_hline(yintercept=bf2_thresh, color="#666666", linetype="dashed") + 
  xlab("log BF1") +
  ylab("log BF2") +
  labs(color = "Result for echolocating lineages:") +
  scale_color_manual(values=corecol(pal="wilke", numcol=2)) +
  guides(color=guide_legend(override.aes=list(alpha=1, size=3))) +
  bartheme() +
  theme(legend.title = element_text(),
        legend.position="bottom")

print(p)
ggsave("fig3a-1.svg", p, width=6, height=5, units="in", dpi=320)
  

# bfs = -20:20
# bf_counts = data.frame("bf1"=c(), "bf2"=c(), "num.acc"=c())
# 
# for(bf1 in bfs){
#   for(bf2 in bfs){
#     cur_acc = nrow(results_raw %>% filter(logbf1 > bf1 & logbf2 > bf2))
#     bf_counts = rbind(bf_counts, data.frame("bf1"=bf1, "bf2"=bf2, "num.acc"=cur_acc))
#     
#   }
# }
# 
# p2 = ggplot(bf_counts, aes(x=bf1, y=bf2, fill=num.acc)) +
#   geom_tile(color="#999999") +
#   scale_x_continuous(expand=c(0,0)) +
#   scale_y_continuous(expand=c(0,0)) +
#   xlab("log BF1 cutoff") +
#   ylab("log BF2 cutoff") +
#   labs(fill="# loci accelerated\nin echolocating lineages") +
#   #scale_fill_gradient(low = corecol(pal="wilke", numcol=1, offset=1), high = corecol(pal="wilke", numcol=1)) +
#   scale_fill_viridis(option="H") +
#   bartheme() +
#   theme(legend.title = element_text(),
#         legend.title.align = 0,
#         legend.position="bottom")
# print(p2)
# ggsave("fig3b.svg", p2, width=6, height=5, units="in", dpi=320)
