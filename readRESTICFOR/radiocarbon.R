####0. Load libraries####

library(lme4)
library(ellipse)
library(Rmisc)
library(lattice)
library(car)
library(carData)
library(ggdendro)
library(AICcmodavg)
library(corrplot)
library(Hmisc)
library(MuMIn)
library(performance)
library(emmeans)
library(sjPlot)
library(ggeffects)
library(lmerTest)
library(ggsignif)
library(broom.mixed)
library(patchwork)
library(tidyverse)

####1. Read and clean radiocarbon data####

bulkD14 <- read.csv("dataRESTICFOR/Delta14C_bulk_soil.csv") %>% 
  separate(id, into = c("crap", "type", "point", "depth"), sep = "_") %>% 
  mutate(site = str_extract(crap, "^[A-Za-z]+")) %>% 
  mutate(crap1 = str_extract(crap, "\\d.*")) %>% 
  mutate(plot = paste0(site, "_", crap1, "_", type)) %>% 
  mutate(unique = paste0(plot, "_", point, "_", depth)) %>% 
  select(-c(F14C, err, err_permil, remark, crap, crap1))

incubations <- read.csv("dataRESTICFOR/Delta14C_incubations.csv") %>% 
  separate(id, into = c("crap", "type"), sep = "_") %>% 
  mutate(site = str_extract(crap, "^[A-Za-z]+")) %>% 
  mutate(crap1 = str_extract(crap, "\\d.*")) %>% 
  mutate(plot = paste0(site, "_", crap1, "_", type)) %>% 
  select(-c(F14C, err, err_permil, crap, crap1))


####2. Exploratory analyses####

#quick dirty analysis for just the upper depth
m_Delta14bulk_sub <- lme4::lmer(Delta14C ~ site * type + (1|site:plot),
                           data = subset(bulkD14, Delta14C >= -200))

summary(m_Delta14bulk_sub)
anova(m_Delta14bulk_sub)
car::Anova(m_Delta14bulk_sub)
r.squaredGLMM(m_Delta14bulk_sub)
windows(12,8)
check_model(m_Delta14bulk_sub)
emmeans(m_Delta14bulk_sub, pairwise ~ site)
plot_model(m_Delta14bulk_sub, type="pred", terms=c("site", "type"))

m_incubations <- lm(Delta14C ~ type *site, data = incubations)

summary(m_incubations)
anova(m_incubations)


####3. Graphs####
bulkD14$type <- recode_factor(bulkD14$type, LONG = "Long", PAST = "Recent")
d15N$site <- recode_factor(bulkD14$site, CT = "Quercus ilex",
                           SET = "Pinus uncinata", SF = "Fagus sylvatica")
bulkD14$site <- factor(bulkD14N$site, levels = c("Quercus ilex", "Fagus sylvatica", 
                                          "Pinus uncinata"))

Delta14C <- ggplot(subset(bulkD14, Delta14C >= -99),
                   aes(x = type, y = Delta14C, fill = type)) +
  geom_boxplot()+
  facet_wrap(~site)+
  scale_fill_manual(values = c("#698B69", "#FFDA7D"))+
  xlab(" ") +
  #ylab(expression(delta^15 * "N (\u2030)")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        text = element_text(size = 21),
        #strip.text = element_blank(),
        axis.line = element_line(color="grey30"),
        axis.ticks.y = element_line(color="black"),
        #axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(vjust = 0.5, size=20),
        axis.title.y = element_blank(),
        plot.margin = margin(1,1,0,0, "cm"))

incubations$type <- recode_factor(incubations$type, LONG = "Long", PAST = "Recent")
incubations$site <- recode_factor(incubations$site, CT = "Quercus ilex",
                           SET = "Pinus uncinata", SF = "Fagus sylvatica")
incubations$site <- factor(incubations$site, levels = c("Quercus ilex", "Fagus sylvatica", 
                                                 "Pinus uncinata"))

windows(12, 6)

incubationsGraph <- ggplot(incubations,
       aes(x = type, y = Delta14C, fill = type)) +
  geom_boxplot() +
  facet_wrap(~site) +
  scale_fill_manual(values = c("#698B69", "#FFDA7D")) +
  xlab(" ") +
  ylab(expression(Delta^14 * "C (\u2030)")) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    
    # PANEL (each facet box)
    panel.background = element_rect(fill = "white", colour = "black"),
    
    # OUTSIDE BACKGROUND (transparent)
    plot.background = element_rect(fill = "transparent", colour = NA),
    
    # Optional: space between panels
    panel.spacing = unit(1, "lines"),
    
    text = element_text(size = 25),
    strip.text = element_blank(),
    
    axis.line = element_line(color = "grey30"),
    axis.ticks.y = element_line(color = "black"),
    axis.title.x = element_blank()
  )

ggsave("incubations.png",
bg = "transparent",width = 15, height = 7, units = "in", dpi = 300)
