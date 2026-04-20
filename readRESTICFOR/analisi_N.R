# 1. Análisi datos de nitrogeno ------------------------------------------------


# 1.1) Datos -------------------------------------------------------------------

#1.1.1 Crga de datos 
dades_N<-readxl::read_excel("dataRESTICFOR/N_data.xlsx")

View(dades_N)

#1.1.1 Limpieza de datos 

#dades_N$Ntotal <- as.numeric(dades_N$Ntotal)

# 1. 2) Representación gráfica  -------------------------------------------------

##### 1.2.1  Gráfico de densidades 

plot(density(na.omit(dades_N$Ntotal)),
     main = "Densitat de Nitrogeno ",
     xlab = "(Ntotal)")

##### 1.2.2  Histograma  

hist(dades_N$Ntotal,
     main = "Histograma de Ntotal",
     xlab = "Ntotal",
     col = "lightblue")

hist(log(dades_N$Ntotal),
     main = "Histograma de log(Ntotal)",
     xlab = "Ntotal",
     col = "lightblue")


##### 1.2.3  QQplot (ver normalidad)

qqnorm(log(dades_N$Ntotal),
       main = "QQ-plot de Nitrogeno ")

qqline(log(dades_N$Ntotal),
       col = "red")


# 2. Modelo lineal mixto  -------------------------------------------------

#2.1 Cargar paquetes necessarios 
library(nlme)

#2.2 Modelos perse: 

m_N_Blanca <- lme4::lmer(log(Ntotal) ~ site * type * SN2 + (1|site:plot),
                         data = dades_N)
summary(m_N_Blanca)
anova(m_N_Blanca)
car::Anova(m_N_Blanca)
r.squaredGLMM(m_N_Blanca)
windows(12, 8)
check_model(m_N_Blanca)
emmeans(m_N_Blanca, pairwise ~ site*type)
emmeans(m_N_sai, pairwise ~ site)
plot_model(m_N_Blanca, type="pred", terms=c("site", "type", "SN2"))

#modelo 1 
modelo1N <- lme(log(Ntotal) ~ site*type*SN2,
               random=~1|site/plot, 
               data=dades_N, 
               na.action = na.omit) 

summary(modelo1N)

anova(modelo1N)

#graficamos el modelo 
library(emmeans)

emm <- emmeans(modelo1N, ~ site * type | SN2)
pairs(emm)
plot(emm,
     by = "SN2",
     xlab = "media estimada Ctotal")


## Extraer residuos del modelo
mis_residuos <- resid(modelo1N)

# Buscar cuál es el valor máximo (que probablemente sea el punto 73)
which.max(mis_residuos)

# Ver la fila completa que tiene ese residuo máximo
dades_N[which.max(mis_residuos), ]


#analisi residuos modelo
library(car)


qqPlot(resid(modelo1N))
plot(resid(modelo2N)~fitted(modelo1N))



#modelo 2 

library(lme4)

modelo2N <- lmer(log(Ntotal) ~ site*type*SN2 +
                  (1 | site/plot),
                data = dades_N, 
                na.action = na.omit)

summary(modelo2N)

anova(modelo2N)

library(car)
Anova(modelo2N, type = "III")  #así nos da el p valor 


# Extraer residuos del modelo
mis_residuos <- resid(modelo2N)

# Buscar cuál es el valor máximo (que probablemente sea el punto 73)
which.max(mis_residuos)

# Ver la fila completa que tiene ese residuo máximo
dades_N[which.max(mis_residuos), ]

emm <- emmeans(modelo2N, ~ site * type | SN2)
pairs(emm)
plot(emm,
     by = "SN2",
     xlab = "media estimada Ntotal")

#ver outliers 

library(car)
outlierTest(modelo2N)


#analisi residuos modelo con qq plot normalidad y con plot podemos ver homostaceidad 
library(car)


qqPlot(resid(modelo2N))
plot(resid(modelo2N)~fitted(modelo2N))
