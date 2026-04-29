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

####1. Read and clean radiocarbon data

bulkD14 <- read.csv("dataRESTICFOR/Delta14C_bulk_soil.csv") %>% 
  separate(id, into = c("crap", "type", "point", "depth"), sep = "_") %>% 
  mutate(site = str_extract(crap, "^[A-Za-z]+")) %>% 
  mutate(crap1 = str_extract(crap, "\\d.*")) %>% 
  mutate(plot = paste0(site, "_", crap1, "_", type)) %>% 
  mutate(unique = paste0(plot, "_", point, "_", depth)) %>% 
  select(-c(F14C, err, err_permil, remark, crap, crap1))

#quick dirty analysis for just the upper depth
m_Delta14bulk_sub <- lme4::lmer(Delta14C ~ site * type + (1|site:plot),
                           data = subset(bulkD14, depth == "0-10"))

summary(m_Delta14bulk_sub)
anova(m_Delta14bulk_sub)
car::Anova(m_Delta14bulk_sub)
r.squaredGLMM(m_Delta14bulk_sub)
windows(12,8)
check_model(m_Delta14bulk_sub)
emmeans(m_Delta14bulk_sub, pairwise ~ site)
eplot_model(m_Delta14bulk_sub, type="pred", terms=c("site", "type"))
