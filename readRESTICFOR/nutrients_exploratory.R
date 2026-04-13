####1. read the data from Jena####

library(tidyverse)
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

jena <- read.csv("dataRESTICFOR/Ctotal_Ntotal_Jena.csv") %>% 
  select(c(Sample.name, SN2, Parameter, Result)) %>%
  mutate(id = paste0(Sample.name, "_", SN2)) %>%
  pivot_wider(names_from = Parameter, values_from = c(Result)) %>%
  mutate(CN_ratio = Ctotal/Ntotal) %>% 
  separate(Sample.name, into = c("plot", "type", "point"), sep = "_") %>% 
  separate(plot, into = c("site", "remove"), sep = "(?=\\d)",
           extra = "merge") %>%
  mutate(plot = paste0(site, "_", type, "_", remove)) %>% 
  select(-c(remove)) %>%
  dplyr::rename(depth = SN2) %>%
  relocate(plot) %>% 
  relocate(c(site, type))

####2. exploratory analyses for total N and C from Jena####
hist(jena$Ctotal)
hist(log(jena$Ctotal))
hist(jena$Ntotal)
hist(log(jena$Ntotal))
hist(jena$CN_ratio)

# log transform C and N content

jena <- jena %>% 
  mutate(lg_Ctotal = log(Ctotal)) %>% 
  mutate(lg_Ntotal = log(Ntotal))

m_Ctotal <- lme4::lmer(log(Ctotal) ~ site * type * depth + (1|site:plot),
                       data = jena)
#always run this before using Anova():
options(contrasts = c("contr.helmert", "contr.poly"))

summary(m_Ctotal)
anova(m_Ctotal)
car::Anova(m_Ctotal)
r.squaredGLMM(m_Ctotal)
check_model(m_Ctotal)
emmeans(m_Ctotal, pairwise ~ site*type)
emmeans(m_Ctotal, pairwise ~ site)
plot_model(m_Ctotal, type="pred", terms=c("site"))
plot_model(m_Ctotal, type="pred", terms=c("type"))
plot_model(m_Ctotal, type="pred", terms=c("site", "type"))
plot_model(m_Ctotal, type="pred", terms=c("site", "type", "depth"))

#check overdispersion:
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m_Ctotal) #there is overdispersion when: ratio > 1

m_Ntotal <- lme4::lmer(log(Ntotal*10) ~ site * type * depth + (1|site:plot),
                       data = jena)

summary(m_Ntotal)
anova(m_Ntotal)
car::Anova(m_Ntotal)
r.squaredGLMM(m_Ntotal)
check_model(m_Ntotal)
emmeans(m_Ntotal, pairwise ~ site*type)
emmeans(m_Ntotal, pairwise ~ site)
plot_model(m_Ntotal, type="pred", terms=c("site"))
plot_model(m_Ntotal, type="pred", terms=c("site", "type", "depth"))
overdisp_fun(m_Ntotal)

m_ratio <- lme4::lmer(CN_ratio ~ site * type * depth + (1|site:plot),
                      data = jena)

summary(m_ratio)
anova(m_ratio)
car::Anova(m_ratio)
r.squaredGLMM(m_ratio)
check_model(m_ratio)
plot_model(m_ratio, type="pred", terms=c("site", "type", "depth"))

####3. Read and pre-process d15N data####

d15N_pre <- read.csv("dataRESTICFOR/d15N_soils.csv") %>% 
  group_by(Muestra) %>% 
  summarise(Nperc = mean(N_perc, na.rm = T),
            d15N_permil = mean(d15N_permil_air, na.rm = T))

d15N <- read.csv("dataRESTICFOR/d15N_soils.csv") %>% 
  separate(Muestra, into = c("site", "crap", "crap2"), sep = "-") %>%
  separate(crap2, into = c("type", "point"), sep = "_") %>% 
  mutate(plot = paste0(site, "_", type, "_", crap)) %>% 
  select(-c(crap, weight_mg, obs)) %>%
  relocate(plot) %>% 
  relocate(c(site, type))

####4. exploratory analyses for total N and d15N from U. Coruña####
hist(d15N$N_perc)
hist(d15N$d15N_permil_air)

N_conc_plot <- 
  ggplot(d15N, aes(x = site, y = N_perc, fill = type)) +
  geom_boxplot() +
  xlab("") +
  ylab("[N] (%)")

m_N_sai <- lme4::lmer(N_perc ~ site * type + (1|site:plot), data = d15N)

summary(m_N_sai)
anova(m_N_sai)
car::Anova(m_N_sai)
r.squaredGLMM(m_N_sai)
check_model(m_N_sai)
emmeans(m_N_sai, pairwise ~ site*type)
emmeans(m_N_sai, pairwise ~ site)
emmeans(m_N_sai, pairwise ~ type)
plot_model(m_N_sai, type="pred", terms=c("site"))
plot_model(m_N_sai, type="pred", terms=c("type"))
plot_model(m_N_sai, type="pred", terms=c("site", "type"))

#check overdispersion:
overdisp_fun(m_N_sai) #there is overdispersion when: ratio > 1

d15N_plot <- 
  ggplot(d15N, aes(x = site, y = d15N_permil_air, fill = type)) +
  geom_boxplot() +
  xlab("") +
  ylab(expression(delta^15 * "N (\u2030)"))

# there are some outliers:
# Above 4 permil: SET-6B-PAST_P3 & SET-9-LONG-P4
# below -3 permil: SET-18-LONG-p3 (only 4 samples in this plot!)

m_d15N <- lme4::lmer(d15N_permil_air ~ site * type + (1|site:plot), data = d15N)

summary(m_d15N)
anova(m_d15N)
car::Anova(m_d15N)
r.squaredGLMM(m_d15N)
check_model(m_d15N)
emmeans(m_d15N, pairwise ~ site*type)
emmeans(m_d15N, pairwise ~ site)
plot_model(m_d15N, type="pred", terms=c("site"))
plot_model(m_d15N, type="pred", terms=c("type"))
plot_model(m_d15N, type="pred", terms=c("site", "type"))
