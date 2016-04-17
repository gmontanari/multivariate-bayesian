library(hglm)
library(sqldf)
library(ggplot2)
data <- read.csv("E:/Users/B146322/Desktop/Fun_estadistica/EST4611_AbortoCrimen.csv")

#REGRESION LINEAL FRECUANTISTA
loghom <- data$lpc_murd
#loghom <- as.numeric(na.omit(loghom))
View(data)
hist(loghom)
N <- length(loghom)
#vemos aqui por ejemplo que la media esta al rededor de -2.5, se ve en el histrograma
#aqui suponemos que la variable se distribuye normal 
#hay info del 72 al 97
sqldf('SELECT 
      AVG(xxprison), AVG(xxpolice), AVG(xxunemp),
      AVG(xxincome), AVG(xxpover),
      STDEV(xxprison), STDEV(xxpolice), 
      STDEV(xxunemp), STDEV(xxincome), 
      STDEV(xxpover)
      FROM data WHERE xxprison IS NOT NULL
      AND xxpolice IS NOT NULL
      AND xxunemp IS NOT NULL
      AND xxincome IS NOT NULL
      AND xxpover IS NOT NULL
      AND year >= 85'
      )

#ESTADISTICOS DE EXPLICATIVAS
avg.year <- sqldf('SELECT year, AVG(lpc_murd) as avg_log_mur 
                  FROM data WHERE lpc_murd IS NOT NULL GROUP BY year')
tab <- data.frame(avg.year$year, exp(avg.year$avg_log_mur))
colnames(tab) <- c('Anio', 'Indice de Homicidios')
p1 <-ggplot(tab, aes(Anio,`Indice de Homicidios`)) + geom_line() 
p1 + ggtitle("Indicide de Homicidios Historico")

#SIN AGRUPAR, DIVIDIDO EN AÑOS
logmurder <- sqldf('SELECT lpc_murd FROM data WHERE lpc_murd IS NOT NULL AND year <85')
logmurder <- data.frame(as.numeric(logmurder$lpc_murd))
colnames(logmurder) <- c('LogHomicidios')
p0 <-ggplot(logmurder, aes(x=LogHomicidios)) + geom_histogram(aes(y=..density..),binwidth=.3, colour="black", fill="white")
p0 + geom_density(alpha=.1, fill="red") + ggtitle("Distribucion de la Media de LogHomicidios")

# Estimador de la media de todos los a??os
data.mean <- mean(logmurder$LogHomicidios) 
data.sd <- sd(logmurder$LogHomicidios)
#estos ya son dos estadisticos

#A PARTIR DEL A??O 81 COMIENZA A DISMINUIR
# NOTE: THESE ESTIMATES ONLY WORK I WE ASSUME WE HAVE A NORMAL DISTRIBUTION!!!!!!!!
alpha <- .05
intervals <- c(data.mean - data.sd*qt(df=N-1, alpha/2, lower.tail=FALSE)/sqrt(N),
               data.mean + data.sd*qt(df=N-1, alpha/2, lower.tail=FALSE)/sqrt(N))

#En todos los casos trabaja con la tasa anual de crecimiento de la incidencia criminal.
exp(intervals)

#A partir de 1991 las tasas de criminalidad, en diferentes dimensiones, 
#experimentaron un decremento sustancial respecto a lo observado en a??os anteriores. 
#Los autores de ese estudio argumentan que la legalizaci??n del aborto a nivel nacional 
#ocurrida en 1973 explica el decremento observado en las tasas de criminalidad 
#observadas 25 a??os despu??s y posteriormente

#EXTRAE Y AGRUPA POR A??OS
data.year<-sqldf('SELECT year, 
              AVG(lpc_murd) as log_mur,
              AVG(xxprison) as prision,
              AVG(xxpolice) as police,
              AVG(xxunemp) as unemp,
              AVG(xxincome) as income,
              AVG(xxpover) as pover
              FROM data 
              WHERE lpc_murd IS NOT NULL
              AND xxprison IS NOT NULL
              AND xxpolice IS NOT NULL
              AND xxunemp IS NOT NULL
              AND xxincome IS NOT NULL
              AND xxpover IS NOT NULL
              GROUP BY year')

data.year <- data.frame(as.numeric(data.year$year), as.numeric(data.year$log_mur),
                        as.numeric(data.year$prision), as.numeric(data.year$police),
                        as.numeric(data.year$unemp), as.numeric(data.year$income), 
                        as.numeric(data.year$pover))

colnames(data.year) <- c('year','log_mur','prision','police','unemp','income','pover')
plot(data.year)

#REGRESION TODAS LAS VARIABLES
cor(data.year,method="pearson")
model<-lm(log_mur ~ prision + police + unemp + income + pover, data.year) 
summary(model)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
#REGRESION VARIABLES SIGNIFICATIVAS
model<-lm(log_mur ~  income + income*unemp + income*pover, data.year) 
summary(model)

model<-lm(log_mur ~  income + income*pover, data.year) 
summary(model)
predict(model)
plot(data.year$log_mur,predict(model))

#EXTRAE SIN NULOS Y SIN AGRUPAR
data1<-sqldf('SELECT year, statenum as state,
              lpc_murd as log_mur,
              xxprison as prision,
              xxpolice as police,
              xxunemp as unemp,
              xxincome as income,
              xxpover as pover
              FROM data 
              WHERE lpc_murd IS NOT NULL
              AND xxprison IS NOT NULL
              AND xxpolice IS NOT NULL
              AND xxunemp IS NOT NULL
              AND xxpover IS NOT NULL')

colnames(data1) <- c('year','state','log_mur','prision','police','unemp','income','pover')
View(data1)
plot(data1)
cor(data1,method="pearson")

model<-lm(log_mur ~ prision + police + unemp + income + pover, data.year) 
summary(model)

#LAS MAS SIGNIFICATIVAS
model<-lm(log_mur ~  income + income*pover, data.year) 
summary(model)

model<-lm(log_mur ~  prision*police + prision, data.year) 
summary(model)


##QUE PASA CON LOS ESTADOS
avg.state <- sqldf('SELECT statenum, log_mur,
                   CASE WHEN year < 85 THEN "bef_85" ELSE "aft_85" END period
                   FROM  
                   (SELECT statenum, year, AVG(lpc_murd) as log_mur    
                   FROM data WHERE lpc_murd IS NOT NULL AND statenum <> 9
                   GROUP BY statenum, year)
                   ')
tab.state <- data.frame(Estado = factor(c(avg.state$statenum)), 
                        Periodo = factor(c(avg.state$period)), 
                        Homicidios = c(exp(avg.state$log_mur)) )
p2 <-ggplot(tab.state, aes(Estado, Homicidios, fill=Periodo)) + geom_bar(stat="identity", position=position_dodge()) 
p2 + ggtitle("Indice de Homicidios por Estado") +  scale_fill_manual(values=c("#666666", "#FF9999"))




##HIERARQUICAL
hier <- hglm(fixed = exp(log_mur) ~ income + income*pover, 
             random = ~1 | state, 
             data = data1, 
             family = Gamma(link = log))
print(summary(hier), print.ranef = TRUE)
plot(hier)

hier2 <- hglm(fixed = log_mur ~ income + income*pover, 
             random = ~1 | state, 
             data = data1, 
             family = gaussian(link = identity))
print(summary(hier2), print.ranef = TRUE)
plot(hier2)

hier3 <- hglm(fixed = exp(log_mur) ~ income + income*pover,
     random = ~ 1|state,
     family = Gamma(link = log),
     disp = ~ prision + police, data = data1)
print(summary(hier3), print.ranef = TRUE)
plot(hier3)

hier4 <- hglm(fixed = log_mur ~ income + income*pover,
              random = ~ 1|state,
              family = gaussian(link = identity),
              disp = ~ prision + police, data = data1)
print(summary(hier4), print.ranef = TRUE)
plot(hier4)


