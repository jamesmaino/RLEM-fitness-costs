library(tidyverse)
library(lme4)

controls = "data/RLEM database for M7, M8 fitness and RMS UM00057.xlsx" %>%
  readxl::read_xlsx(sheet="data") %>% 
  mutate(trtcode = ifelse(is.na(trtcode), 1, trtcode)) %>%
  filter(trtcode == 1) %>% 
  distinct(plot) %>%
  pull(plot)

d = readxl::read_xlsx("data/Ant_RLEM_MiSeq_v3_svs_JM.xlsx", sheet="TRANSPOSED") %>% 
  rename(Site = "Field Site") %>% 
  mutate(Date = as.Date(Date, format="%d/%m/%Y")) %>%
  filter(Plot %in% controls)
code = readxl::read_xlsx("data/Ant_RLEM_MiSeq_v3_svs_JM.xlsx", sheet="otu_to_R") %>% 
  rename(name = "otu")


dsum = d %>% 
  select(-S, -G, -TTC) %>%
  pivot_longer(R:OTU_19) %>% 
  group_by(Site, Date, Plot, name) %>% 
  summarise(value = mean(value), .groups="drop") %>%
  group_by(Site, Date, name) %>% 
  summarise(
    N=n(),
    mean_val = mean(value),
    se_val = sd(value)/sqrt(n()),
    .groups="drop")

# plot OTUs

dsum %>% 
  filter(grepl("OTU", name)) %>%
  mutate(name = ifelse(name %in% paste0("OTU_", 1:5), name, "other")) %>%
  group_by(Site, Date, name) %>% 
  summarise(mean_val = sum(mean_val), .groups="drop") %>%
  mutate(name = reorder(name, as.numeric(str_extract(name, "\\d+")))) %>%
  ggplot(aes(Date, mean_val, fill=name)) +
  geom_area(colour="black") + 
  facet_wrap(~Site) + 
  theme_bw() +
  coord_cartesian(expand = FALSE)
ggsave("plots/explore_amplicon/OTUs.png", height=4, width=10)



# plot resistance genotype for arthur river
d1 = d %>% 
  select(-S, -G, -TTC) %>%
  pivot_longer(R:OTU_19) %>% 
  left_join(code) %>%
  drop_na(R) %>% 
  # filter(Date == "2017-07-31") %>%
  filter(Barcode != "plate2b06_S138") %>% # plate2b06_S138 overlaps with plate7e02_S525
  filter(Barcode != "plate1b01_S2") %>% # plate6b01_S482 overlaps with plate1b01_S2
  filter(Barcode != "plate1f01_S6") %>% # plate6f01_S486 overlaps with plate1f01_S6
  filter(Barcode != "plate1e01_S5") %>% # plate6e01_S485 overlaps with plate1e01_S5
  filter(Barcode != "plate1c01_S3") %>% # plate6c01_S483 overlaps with plate1c01_S3 
  filter(Barcode != "plate1d01_S4") %>% # plate6d01_S484 overlaps with plate1d01_S4 
  filter(Site == "Arthur River") %>%  
  group_by(Site, Date, Plot, Replicate, R) %>% 
  summarise(value = sum(value), .groups="drop") %>% 
  mutate(time_d = Date - min(Date)) %>%
  filter(R != "TTG")
  
d1 %>% 
  group_by(Site, Date, Plot, R) %>% 
    summarise(
    N=n(),
    mean_val = mean(value),
    .groups="drop") %>%
   group_by(Site, Date, R) %>% 
    summarise(
    N=n(),
    se_val = sd(mean_val)/sqrt(N),
    mean_val = mean(mean_val),
    .groups="drop") %>%
  ggplot(aes(Date, mean_val, fill=R, color=R)) +
  # geom_area() +
  # geom_point() +
    geom_line() +
  geom_ribbon(color=FALSE,aes(ymin=mean_val-se_val, ymax=mean_val+se_val), alpha=0.4) +
  theme_bw() + 
  ylab("Mean allele frequency") +
  coord_cartesian(expand = FALSE) + 
  theme(axis.text.x = element_text(angle=45, hjust=1))
ggsave("plots/explore_amplicon/Rgenotypes.png", height=4, width=6)


# check whether there is support to separate TTT and TTC
m1 = lm((value) ~ R/time_d + Plot - 1, data=d1)
m1 = lmer((value) ~ R*time_d + (1|Plot), data=filter(d1, value!=0))
summary(m1)
car::Anova(m1)

# plot SNPs
dsum %>% 
  filter(!grepl("OTU", name)) %>%
  ggplot(aes(Date, mean_val, color=name, fill=name)) +
  geom_line() +
  geom_ribbon(color=FALSE,aes(ymin=mean_val-se_val, ymax=mean_val+se_val), alpha=0.4) +
  geom_point() +
  facet_wrap(~Site) + 
  theme_bw() 
ggsave("plots/explore_amplicon/SNPs.png", height=4, width=10)
