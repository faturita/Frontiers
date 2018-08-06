##fijo el directorio de trabajo
setwd("C:/Users/Debie/Dropbox/rodrigo")
getwd()
####cargo las librerías
library(nlme)
library(car)
library(ggplot2)
library(grid)
library(datasets)
library(RColorBrewer)
library(extrafont)
library(faraway)
library(effects)
### defino las variables de la base
method=factor(c(rep("M1",16),rep("M2",16),rep("M3",16)))
ind=(c(rep(seq(1:16),3)))
group=(1*(ind<9)+2*(ind>8))
Ind=factor(ind)
Group=factor(group)
percent=c(0.35,0.85,0.25,0.55,0.4,0.6,0.8,0.95,0.4,0.3,0.4,0.45,0.6,0.5,0.7,0.5,
       0.45,0.3,0.65,0.4,0.35,0.35,0.6,0.9,0.65,0.15,0.5,0.4,0.3,0.35,0.25,0.35,
       0.4,0.5,0.55,0.5,0.45,0.7,0.35,0.95,0.4,0.1,0.25,0.2,0.2,0.3,0.3,0.2)
##le doy formato de dataframe
rodrigo=data.frame(method,percent,group,ind,Group,Ind)


###resumenes por grupo y por método
tapply(percent, method, summary)
tapply(percent, group, summary)

####visualizamos la información
p1 <- ggplot(rodrigo, aes(x = Group, y = percent,fill=method)) +
  geom_boxplot(   alpha = 0.77) +
  scale_y_continuous(name = "% Success", breaks = seq(0, 1, 0.25), limits=c(0, 1.1)) +
  scale_x_discrete(name = "Group") +
  ggtitle("Success by group and Method") +
  theme_bw()  +
  scale_fill_brewer(palette = "Accent")+  facet_grid(. ~method)
p1 
ggsave(file="C:/Users/Debie/Dropbox/rodrigo/plot1.pdf", p1)


####modelo sin efectos aleatorios de individuo

###las interacciones
par(mfrow=c(1,2))
interaction.plot( method, Group, percent, legend = T)
interaction.plot(Group, method, percent, legend = T)
perc.aov2=aov(percent ~ Group * method)
summary(perc.aov2)

## descartamos las interacciones por no ser estadisticamente significativas
rodri.eff=effect("Group*method",model.rod,confidence.level=0.90)
plot(rodri.eff,main="")
perc.aov1=aov(percent ~ Group + method)
summary(perc.aov1)
model.tables(perc.aov1, type = "means") # visualizamos el efecto por método y por grupo
plot(perc.aov1)
###Comparaciones a posteriori
TukeyHSD(perc.aov1, "Group",ordered=T)
TukeyHSD(perc.aov1, "method",ordered=T)
p3=plot(TukeyHSD(perc.aov1, "method",ordered=T))
ggsave(file="C:/Users/Debie/Dropbox/rodrigo/plot3.pdf",plot(TukeyHSD(perc.aov1, "method",ordered=T)))
###modelo de análisis de la varianza con efectos aleatorios 
lme.rodri1 <- lme(percent ~ Group + method, random=~1|Ind, data=rodrigo)
anova(lme.rodri1)
summary(lme.rodri1)
plot(lme.rodri1)

#análisis diagnóstico del modelo
# normalidad de los residuos
shapiro.test(residuals(lme.rodri1))

qqPlot(residuals(lme.rodri1),ylab = "residuals",
       col = "coral",pch = 19,
       col.lines = "cadetblue")

ggsave(file="C:/Users/Debie/Dropbox/rodrigo/plot2.pdf",
       qqPlot(residuals(lme.rodri1),ylab = "residuals",
              col = "coral",pch = 19,
              col.lines = "cadetblue"))
#homocedasticidad de los residuos
leveneTest(percent , Group )
leveneTest(percent , method )
