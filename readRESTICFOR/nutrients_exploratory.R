####1. read the data####

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

####2. exploratory analyses####
hist(jena$Ctotal)
hist(log(jena$Ctotal))
hist(jena$Ntotal)
hist(log(jena$Ntotal))
hist(jena$CN_ratio)

# log transform C and N content

jena <- jena %>% 
  mutate(lg_Ctotal = log(Ctotal)) %>% 
  mutate(lg_Ntotal = log(Ntotal))

m_Ctotal <- lme4::lmer(lg_Ctotal ~ site * type * depth + (1|site:plot),
                       data = jena)
#always run this before using Anova():
options(contrasts = c("contr.helmert", "contr.poly"))

summary(m_Ctotal)
anova(m_Ctotal)
car::Anova(m_Ctotal)
r.squaredGLMM(m_Ctotal)
check_model(m_Ctotal)
emmeans(m_Ctotal, pairwise ~ site*type)
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

m_Ntotal <- lme4::lmer(lg_Ntotal ~ site * type * depth + (1|site:plot),
                       data = jena)

summary(m_Ntotal)
anova(m_Ntotal)
car::Anova(m_Ntotal)
r.squaredGLMM(m_Ntotal)
check_model(m_Ntotal)
emmeans(m_Ntotal, pairwise ~ site*type)
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
