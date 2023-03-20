library(tidyverse)
library(lubridate)
library(cowplot) 

datafile = "data/UM00057_field fitness 2017 SA only_Aug 2020.xlsx"

codes = datafile %>% 
  readxl::read_xlsx(sheet = 1) %>% 
  select(trtcode, Year, Treatment)

d = datafile %>% 
  readxl::read_xlsx(sheet = 2) %>% 
  select(date, location, month, plot, trtcode, 
         `abundance per m2`, kdr_percent.mean) %>% 
  mutate(Year = year(date)) %>% 
  left_join(codes, by = c("trtcode", "Year"))

dsum = d %>% 
  filter(trtcode %in% c(1,2,4)) %>%
  group_by(month, Treatment) %>%
  summarise(
    abund    = mean(`abundance per m2`),
    abund_se = sd(`abundance per m2`)/sqrt(n()), 
    kdr    = mean(kdr_percent.mean),
    kdr_se = sd(kdr_percent.mean)/sqrt(n())
  ) %>% 
  mutate(month = factor(month, 
                        levels = c("June", "August","September","October"))) %>%
  mutate(Treatment = factor(Treatment, 
                        levels = c("Untreated Control", "Bifenthin x1", "Bifenthin x2")))

  
# abund
p1 = dsum %>% 
  ggplot(aes(Treatment, abund)) + 
  geom_bar(stat = "identity", alpha=0.5) + 
  geom_errorbar(aes(ymin = abund-abund_se, ymax = abund+abund_se),
                width=0.2) +
  facet_grid(~month) + 
  xlab("Treatment") + 
  ylab(expression("Average number of mites per m"^2)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 4000))
p1
ggsave("plots/Fig5a_SPtreatment_abund.jpeg", width=10, height = 4, scale = 0.9)  

# stats
ds = d %>% 
  filter(trtcode %in% c(1,2,4)) %>% 
  mutate(trtcode = factor(trtcode))

# There was an initial decline in mite numbers following the first chemical treatment in June 2017 (compared with the untreated controls), with a significant difference between the controls and each of the treatments (t = ), 

lm(log(`abundance per m2`) ~ trtcode, 
   data=filter(ds, month == "August")) %>% 
  summary()


# although the 1x and 2x treatments were similar (t =). 

lm(log(`abundance per m2`) ~ trtcode, 
   data=filter(ds, month == "August") %>% 
     mutate(trtcode = relevel(trtcode, ref = 2))) %>% 
  summary()


# However, mite numbers rebounded quickly, surpassing densities within the control treatments, although not significantly (t = etc) 

lm(log(`abundance per m2`) ~ trtcode, 
   data=filter(ds, month == "September")) %>% 
  summary()


 # kdr
p2 = dsum %>% 
  ggplot(aes(Treatment, kdr)) + 
  geom_bar(stat = "identity", alpha=0.5) + 
  geom_errorbar(aes(ymin = kdr-kdr_se, ymax = kdr+kdr_se),
                width=0.2) +
  facet_grid(~month) + 
  xlab("Treatment") + 
  ylab(expression("Average "*italic(kdr)*' allele frequency (%)')) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 110))
p2
ggsave("plots/Fig5b_SPtreatment_kdr.jpeg", width=10, height = 4, scale = 0.9)  
