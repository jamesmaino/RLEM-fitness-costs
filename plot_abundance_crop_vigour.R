library(tidyverse)

# archived spreadsheet
d0 = "data/RLEM database for M7, M8 fitness and RMS UM00057.xlsx" %>%
  readxl::read_xlsx(sheet="data") %>% 
  mutate(damage = ifelse(is.na(damage), 0, damage))

d  = d0 %>% 
  # remove treatment plots from Tintinara
  mutate(trtcode = ifelse(is.na(trtcode), 1, trtcode)) %>%
  filter(trtcode == 1) %>% 
  
  # samples were taken over multiple dates so rough date is used which is the last date of sampling for a given field trip
  mutate(rough_date = as.Date(rough_date)) %>%
  
  # get mean across quadrats where multiple samples are available
  dplyr::select(rough_date,	location,	month, plot,	
                trtcode, damage, contains("abund.q")) %>%
  gather(var, val, -rough_date,	-location,	-month,	-plot,	-trtcode,	-damage) %>%    separate(var, c('var', 'quadrat'), sep = '\\.q', convert = TRUE) %>% 
  filter(!is.na(quadrat)) %>%
  filter(!is.na(val)) %>%
  spread(var, val, convert = TRUE) %>%
  group_by(rough_date, location, month, plot, trtcode) %>% 
  summarise(abund=mean(abund), damage=damage[1], .groups="drop") %>%
  # convert 30cm x 30 cm abund to m2
  mutate(abund_m2 = abund/(0.3^2))

# make summary table of abundance and damage
dsum = d %>% 
  filter(rough_date != "2019-10-22") %>%
  mutate(month_year = format(rough_date, "%b %Y")) %>%
  mutate(month_year = fct_reorder(month_year, rough_date)) %>%
  filter(location %in% c("Tintinara", "Arthur River")) %>% 
  group_by(rough_date, month_year, location, month) %>% 
  summarise(
    abund_m2_mu=mean(abund_m2), 
    abund_m2_se=sd(abund_m2)/sqrt(n()), 
    damage_mu=mean(damage), 
    damage_se=sd(damage)/sqrt(n()), 
    .groups="drop") %>% 
  arrange(location)
write_csv(dsum, "plots/summary_table_abund_damage.csv")

# 1. RLEM abundances +- SE (extrapolated to m2) by sampling date for Tintanara (excluding data for 22 October 2019 which was too late in the year)
# 2. RLEM abundances +- SE (extrapolated to m2) by sampling date for Arthur River (excluding data for 22 October 2019 which was too late in the year)

dsum %>%
  mutate(ymin = abund_m2_mu-abund_m2_se, 
         ymax = abund_m2_mu+abund_m2_se) %>%
  ggplot(aes(month_year, abund_m2_mu)) + 
  # geom_point() +
  geom_bar(stat="identity", alpha=0.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.3) +
  facet_wrap(~location, scales = "free_x") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90)) +
  # scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(drop=TRUE) +
  xlab("Date") +
  ylab(expression("Average number of mites per m"^2))
ggsave('plots/Fig1_abundance_mu_se.jpeg', width=6, height=3.5, dpi = 300)


# 3. Mite feeding damage +- SE (score out of 10) by sampling date for Tintanara (excluding data for 22 October 2019 which was too late in the year)
# 4. Mite feeding damage +- SE(score out of 10) by sampling date for Arthur River (excluding data for 22 October 2019 which was too late in the year). NB: whenever Owain entered NA is meant zero damage

dsum %>%
  mutate(ymin = damage_mu-damage_se, ymax = damage_mu+damage_se) %>%
  ggplot(aes(month_year, damage_mu)) + 
  # geom_point() +
  geom_bar(stat="identity", alpha=0.5) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=0.3) +
  facet_wrap(~location, scales = "free_x") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90)) +
  # scale_y_continuous(expand=c(0,0)) +
  ylim(0, 10) +
  xlab("Date") +
  ylab(expression("Average feeding damage score"))
ggsave('plots/Fig2_damage_mu_se.jpeg', width=6, height=3.5, dpi = 300)

