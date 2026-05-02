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

####1. read the data from Jena####

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
#####2.1 analyses of C from Jena#####

m_Ctotal <- lme4::lmer(log(Ctotal) ~ site * type * depth + (1|site:plot),
                       data = jena)
#quick dirty analysis for just the upper depth
m_Ctotal_sub <- lme4::lmer(log(Ctotal) ~ site * type + (1|site:plot),
                       data = subset(jena, depth == "0-10 cm"))
#without log transforming
m_Ctotal_raw <- lme4::lmer(Ctotal ~ site * type * depth + (1|site:plot),
                       data = jena)
#quick dirty analysis for just the upper depth
m_Ctotal_sub_raw <- lme4::lmer(Ctotal ~ site * type + (1|site:plot),
                           data = subset(jena, depth == "0-10 cm"))
#always run this before using Anova():
options(contrasts = c("contr.helmert", "contr.poly"))

summary(m_Ctotal)
anova(m_Ctotal)
car::Anova(m_Ctotal)
r.squaredGLMM(m_Ctotal)
windows(12, 8)
check_model(m_Ctotal)
emmeans(m_Ctotal, pairwise ~ site*type)
emmeans(m_Ctotal, pairwise ~ site*type*depth)
emmeans(m_Ctotal, pairwise ~ site)
plot_model(m_Ctotal, type="pred", terms=c("site"))
plot_model(m_Ctotal_raw, type="pred", terms=c("type"))
plot_model(m_Ctotal_raw, type="pred", terms=c("site", "type"))
plot_model(m_Ctotal_raw, type="pred", terms=c("site", "type", "depth"))

car::Anova(m_Ctotal_sub)
windows(12, 8)
check_model(m_Ctotal_sub)
emmeans(m_Ctotal_sub, pairwise ~ site*type)
emmeans(m_Ctotal_sub, pairwise ~ site)
plot_model(m_Ctotal_sub, type="pred", terms=c("site"))
plot_model(m_Ctotal_sub, type="pred", terms=c("type"))
plot_model(m_Ctotal_sub, type="pred", terms=c("site", "type"))

car::Anova(m_Ctotal_raw)
windows(12, 8)
check_model(m_Ctotal_raw)
emmeans(m_Ctotal_raw, pairwise ~ site*type)
emmeans(m_Ctotal_raw, pairwise ~ site)
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("site"))
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("type"))
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("site", "type", "depth"))


car::Anova(m_Ctotal_sub_raw)
windows(12, 8)
check_model(m_Ctotal_sub_raw)
emmeans(m_Ctotal_sub_raw, pairwise ~ site*type)
emmeans(m_Ctotal_sub_raw, pairwise ~ site)
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("site"))
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("type"))
plot_model(m_Ctotal_sub_raw, type="pred", terms=c("site", "type"))

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

#####2.2 analyses of N from Jena#####
m_Ntotal <- lme4::lmer(log(Ntotal*10) ~ site * type * depth + (1|site:plot),
                       data = jena)
#quick dirty analysis for just the upper depth
m_Ntotal_sub <- lme4::lmer(log(Ntotal) ~ site * type + (1|site:plot),
                           data = subset(jena, depth == "0-10 cm"))

summary(m_Ntotal)
anova(m_Ntotal)
car::Anova(m_Ntotal)
r.squaredGLMM(m_Ntotal)
check_model(m_Ntotal)
emmeans(m_Ntotal, pairwise ~ site*type)
emmeans(m_Ntotal_sub, pairwise ~ site)
plot_model(m_Ntotal, type="pred", terms=c("site"))
plot_model(m_Ntotal_sub, type="pred", terms=c("site", "type"))
overdisp_fun(m_Ntotal)

car::Anova(m_Ntotal_sub)
windows(12, 8)
check_model(m_Ctotal_sub)
emmeans(m_Ntotal_sub, pairwise ~ site*type)
emmeans(m_Ctotal_sub, pairwise ~ site)
plot_model(m_Ntotal_sub, type="pred", terms=c("site"))
plot_model(m_Ntotal_sub, type="pred", terms=c("type"))
plot_model(m_Ntotal_sub, type="pred", terms=c("site", "type"))

#####2.3. analyses of C:N ratio#####
m_ratio <- lme4::lmer(CN_ratio ~ site * type * depth + (1|site:plot),
                      data = jena)

summary(m_ratio)
anova(m_ratio)
car::Anova(m_ratio)
r.squaredGLMM(m_ratio)
check_model(m_ratio)
plot_model(m_ratio, type="pred", terms=c("site", "type", "depth"))

####3. Read and pre-process d15N data####

d15N <- read.csv("dataRESTICFOR/d15N_soils.csv") %>% 
  #some samples were measured twice
  group_by(Muestra) %>% 
  summarize(Nperc = mean(N_perc, na.rm = T),
            d15N_permil = mean(d15N_permil_air, na.rm = T)) %>% 
  separate(Muestra, into = c("site", "crap", "crap2"), sep = "-") %>%
  separate(crap2, into = c("type", "point"), sep = "_") %>% 
  mutate(plot = paste0(site, "_", type, "_", crap)) %>% 
  select(-c(crap)) %>%
  relocate(plot) %>% 
  relocate(c(site, type))

####4. exploratory analyses for total N and d15N from U. Coruña####
hist(d15N$Nperc)
hist(d15N$d15N_permil)

N_conc_plot <- 
  ggplot(d15N, aes(x = site, y = Nperc, fill = type)) +
  geom_boxplot() +
  xlab("") +
  ylab("[N] (%)")

m_N_sai <- lme4::lmer(Nperc ~ site * type + (1|site:plot), data = d15N)

summary(m_N_sai)
anova(m_N_sai)
car::Anova(m_N_sai)
r.squaredGLMM(m_N_sai)
windows(12, 8)
check_model(m_N_sai)
emmeans(m_N_sai, pairwise ~ site*type)
emmeans(m_N_sai, pairwise ~ site)
emmeans(m_N_sai, pairwise ~ type)
plot_model(m_N_sai, type="pred", terms=c("site"))
plot_model(m_N_sai, type="pred", terms=c("type"))
plot_model(m_N_sai, type="pred", terms=c("site", "type"))

d15N_plot <- 
  ggplot(d15N, aes(x = site, y = d15N_permil, fill = type)) +
  geom_boxplot() +
  xlab("") +
  ylab(expression(delta^15 * "N (\u2030)"))

m_d15N <- lme4::lmer(d15N_permil ~ site * type + (1|site:plot), data = d15N)

summary(m_d15N)
anova(m_d15N)
car::Anova(m_d15N)
r.squaredGLMM(m_d15N)
windows(12, 8)
check_model(m_d15N)
emmeans(m_d15N, pairwise ~ site*type)
emmeans(m_d15N, pairwise ~ site)
plot_model(m_d15N, type="pred", terms=c("site"))
plot_model(m_d15N, type="pred", terms=c("type"))
plot_model(m_d15N, type="pred", terms=c("site", "type"))

####5. Method comparison (Jena vs. IRNAS)####

Ctotal <- read.csv("dataRESTICFOR/C_Total_Resticfor_IRNAS.csv")

comp <- Ctotal %>% 
  filter(purpose == "method_comparison") %>% 
  separate(MUESTRA, into = c("plot", "type", "point", "depth"), sep = "_") %>%
  mutate(site = str_extract(plot, "^[A-Za-z]+")) %>% 
  mutate(crap = str_extract(plot, "\\d.*")) %>% 
  mutate(depth = ifelse(depth == "10 cm", "0-10 cm", "10-20 cm")) %>% 
  mutate(id = paste0(plot, "_", type, "_", point, "_", depth)) %>% 
  mutate(plot = paste0(site, "_", type, "_", crap)) %>%
  select(c(id, C_total_perc)) %>% 
  rename(C_total_perc_IRNAS = C_total_perc) %>% 
  left_join(jena, by = "id")
 
summary(lm(Ctotal ~ C_total_perc_IRNAS, data = comp))
# there is an outlier for the Jena data base, try removing
summary(lm(Ctotal ~ C_total_perc_IRNAS, data = subset(comp, Ctotal <= 17)))
confint(lm(Ctotal ~ C_total_perc_IRNAS, data = subset(comp, Ctotal <= 17)))
# R2 improves (0.96), but the slope is significantly different from 1
# that is values from IRNAS are consistently lower than from Jena

C_total_comp <- 
  ggplot(subset(comp, Ctotal <= 17), aes(x = C_total_perc_IRNAS, y = Ctotal)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  xlab("[C] (%) IRNAS") +
  ylab("[C] (%) Jena")

####6. Analyses of total C from soil samples 0-20 (IRNAS)####

####7. Graphs####

jena$site <- recode_factor(jena$site, CT = "Quercus ilex",
                           SET = "Pinus uncinata", SF = "Fagus sylvatica")
jena$site <- factor(jena$site, levels = c("Quercus ilex", "Fagus sylvatica",
                                          "Pinus uncinata"))
jena$type <- recode_factor(jena$type, LONG = "Long", PAST = "Recent")

windows(16, 8)
C_jena_0to10 <- ggplot(subset(jena, depth == "0-10 cm"),
                       aes(x = type, y = Ctotal, fill = type))+
  geom_boxplot() +
  facet_wrap(~site) +
  scale_fill_manual(values = c("#698B69", "#FFDA7D")) +
  xlab(" ") +
  ylab("Soil C (%)") +
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

ggsave("carbon_0to10.png",
       bg = "transparent",width = 15, height = 7, units = "in", dpi = 300)



N_0to10cm_graph <-
  ggplot(subset(jena, depth == "0-10 cm"),
         aes(x = site, y = Ntotal, fill = type)) +
  geom_boxplot(aes(fill = type), width = 0.7) +
  xlab(" ") +
  ylab("Soil [N] (%)") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 21),
        axis.line = element_line(color="grey30"),
        axis.ticks.y = element_line(color="black"),
        axis.ticks.x =element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 0.5, size=20),
        plot.margin = margin(1,0,0,0, "cm")) 

windows(12, 8)
cowplot::plot_grid(C_0to10cm_graph, N_0to10cm_graph)

d15N$type <- recode_factor(d15N$type, LONG = "Long", PAST = "Recent")
d15N$site <- recode_factor(d15N$site, CT = "Quercus ilex",
                          SET = "Pinus uncinata", SF = "Fagus sylvatica")
d15N$site <- factor(d15N$site, levels = c("Quercus ilex", "Fagus sylvatica", 
                                          "Pinus uncinata"))

windows(20, 9)

Ngraph<- ggplot(d15N, aes(x = type, y = Nperc, fill = type)) +
  geom_boxplot() +
  facet_wrap(~site) +
  scale_fill_manual(values = c("#698B69", "#FFDA7D")) +
  xlab(" ") +
  ylab("Soil N (%)") +
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


d15Ngraph <- ggplot(d15N, aes(x = type, y = d15N_permil, fill = type))  +
  geom_boxplot() +
  facet_wrap(~site) +
  scale_fill_manual(values = c("#698B69", "#FFDA7D")) +
  xlab(" ") +
  ylab(expression(delta^15 * "N (\u2030)")) +
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

cowplot::plot_grid(Ngraph, d15Ngraph, nrow =2)
ggsave("Nitrogen_d15N.png",
       bg = "transparent",width = 15, height = 7, units = "in", dpi = 300)
