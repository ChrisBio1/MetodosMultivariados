library(vegan)

especies <- read.csv("Datos/especies.csv", header = T)
indices <- read.csv("Datos/indices.csv", header = T)
ubicacion <- read.csv("Datos/ubicacion.csv", header = T)

especiesA <- especies[c(-1:-3)]
indicesA <- indices[c(5:10)]


### Análisis de conglomerados --------------------------------------
ob1 <- data.frame(a = c(0, 1, 4, 5),
                  b = c(0, 2, 3, 6),
                  c = c(0, 1, 2, 7),
                  d = c(0, 1, 2, 5))
vegdist(ob1, "euclidea")

dev.new()
plot(hclust(vegdist(ob1, "euclidea"), "single"))
abline(h = 2.64, lty = 3)
abline(h = 3.46, col = "red", lty = 3)
abline(h = 6.63, col = "blue", lty = 3)

ob2 <- data.frame(a = rpois(4, 5),
                  b = rpois(4, 7),
                  c = rpois(4, 6),
                  d = rpois(4, 5))

plot(hclust(vegdist(ob2, "euclidea"), "single"))
# Método de aglomeramiento del vecino más cercano

especiesC <- rowsum(especiesA, paste0(especies$punto, especies$habitat))
dim(especiesC)
espN <- decostand(especiesC, "normalize")
espeNEU <- vegdist(espN, "euc")
espeNEU.VC <- hclust(espeNEU, method="single")
# options(width = 200)
?hclust()
plot(espeNEU.VC)
summary(espeNEU.VC)
dev.new()
plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"), "complete"))

# Método de aglomeramiento del vecino más lejano
espeNEU.VL <- hclust(espeNEU, method="complete")
plot(espeNEU.VL)

plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"),"complete"))

# Método de aglomeramiento UPGMA
espeNEU.UPGMA <- hclust(espeNEU, method="average")
plot(espeNEU.UPGMA)

plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"),"average"))

# Método de aglomeramiento por centroide
espeNEU.CEN <- hclust(espeNEU, method="centroid")
plot(espeNEU.CEN)

plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"),"centroid"))

# Método de aglomeramiento por varianza mínima de Ward
espeNEU.WA <- hclust(espeNEU, method="ward.D2")
plot(espeNEU.WA)

plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"),"ward.D2"))

### Paréntesis Autos
autos <- read.csv("Autos/autos.csv", header = T, sep = "\t")
m1 <- head(autos, 8)[c(5:38)]

plot(hclust(vegdist(m1 ,"euclidean"), "ward.D"))
plot(hclust(vegdist(decostand(m1, "norm") ,"euclidean"), "ward.D"))
###

### Correlación cofenética -------------------------------------
# Vecino más cercano
espeCofC <- cophenetic(espeNEU.VC)
cofC <- cor(espeNEU, espeCofC)
cor.test(espeNEU,espeCofC)

# Vecino más lejano
espeCofL <- cophenetic(espeNEU.VL)
cofL <- cor(espeNEU, espeCofL)

# UPGMA
espeCofU <- cophenetic(espeNEU.UPGMA)
cofU <- cor(espeNEU, espeCofU)

# Centroide
espeCofCEN <- cophenetic(espeNEU.CEN)
cofCE <- cor(espeNEU, espeCofCEN)

# Ward
espeCofW <- cophenetic(espeNEU.WA)
cofW <- cor(espeNEU, espeCofW)

data.frame(cofC, cofL, cofU, cofCE, cofW)

### Graficación de correlación de matrices
dev.new()
par(mfrow=c(2,3))
plot(espeNEU, espeCofC,
     xlab="Distancia Cuerda", 
     ylab="Distancia cofenética", 
     asp=1, 
     xlim=c(0,sqrt(1)), ylim=c(0,sqrt(1)),
     main=c("Vecino más cercano", 
            paste("Correlación cofenética =",
                                    round(cofC,3))))
abline(0,1)
lines(lowess(espeNEU, espeCofC), col="red")

plot(espeNEU, espeCofL,
     xlab="Distancia Cuerda", 
     ylab="Distancia cofenética", 
     asp=1, 
     xlim=c(0,sqrt(1)), ylim=c(0,sqrt(1)),
     main=c("Vecino más lejano", 
            paste("Correlación cofenética =",
                                    round(cofL,3))))
abline(0,1)
lines(lowess(espeNEU, espeCofL), col="red")

plot(espeNEU, espeCofU,
     xlab="Distancia Cuerda", 
     ylab="Distancia cofenética", 
     asp=1, 
     xlim=c(0,sqrt(1)), ylim=c(0,sqrt(1)),
     main=c("UPGMA", 
            paste("Correlación cofenética =",
                                    round(cofU,3))))
abline(0,1)
lines(lowess(espeNEU, espeCofU), col="red")


plot(espeNEU, espeCofCEN,
     xlab="Distancia Cuerda", 
     ylab="Distancia cofenética", 
     asp=1, 
     xlim=c(0,sqrt(1)), ylim=c(0,sqrt(1)),
     main=c("Centroide", 
            paste("Correlación cofenética =",
                                    round(cofCE,3))))
abline(0,1)
lines(lowess(espeNEU, espeCofCEN), col="red")

plot(espeNEU, espeCofW,
     xlab="Distancia Cuerda", 
     ylab="Distancia cofenética", 
     asp=1, 
     xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Ward", 
            paste("Correlación cofenética =",
                                    round(cofW,3))))
abline(0,1)
lines(lowess(espeNEU, espeCofW), col="red")

plot(espeNEU, espeNEU,
     xlab="Distancia Cuerda", 
     ylab="Distancia Cuerda", 
     asp=1, 
     xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Cuerda", 
            paste("Correlación nula =",
                                    round(cor(espeNEU, espeNEU),3))))
abline(0,1)
lines(lowess(espeNEU, espeNEU), col="red")


# Distancias de Gower (1983)
sum((espeNEU- cofC)^2) # Vecino más cercano 
sum((espeNEU- cofL)^2) # Vecino más lejano
sum((espeNEU- cofU)^2) # UPGMA
sum((espeNEU- cofCE)^2) # Centroide
sum((espeNEU- cofW)^2) # Ward


### Elección del número de Conglomerados ---------------------------------

# Gráficas de fusión de niveles
plot(hclust(vegdist(decostand(especiesC, "norm"), "euclidea"), "average"))
# Número de conglomerados por vecino más cercano
dev.new()
plot(espeNEU.VC$height, nrow(especiesC):2, type="S", 
     main="Fusion levels - Cuerda - V. cercano", 
     ylab="k (no. conglomerados)", 
     xlab="h (altura del nodo)", 
     col="grey")
text(espeNEU.VC$height, 
     nrow(especiesC):2, 
     nrow(especiesC):2, 
     col="red", cex=0.8)

# Número de conglomerados por vecino más lejano
plot(espeNEU.VL$height, nrow(especiesC):2, type="S", 
     main="Fusion levels - Cuerda - V. lejano", 
     ylab="k (no. conglomerados)", 
     xlab="h (altura del nodo)", 
     col="grey")
text(espeNEU.VL$height, 
     nrow(especiesC):2, 
     nrow(especiesC):2, 
     col="red", cex=0.8)

# Número de conglomerados por UPGMA
plot(espeNEU.UPGMA$height, nrow(especiesC):2, type="S", 
     main="Fusion levels - Cuerda - UPGMA", 
     ylab="k (no. conglomerados)", 
     xlab="h (altura del nodo)", 
     col="grey")
text(espeNEU.UPGMA$height, 
     nrow(especiesC):2, 
     nrow(especiesC):2, 
     col="red", cex=0.8)
plot(hclust(vegdist(decostand(especiesC, "norm"), "euclidea"), "average"))

# Número de conglomerados por centroide
plot(espeNEU.CEN$height, nrow(especiesC):2, type="S", 
     main="Fusion levels - Cuerda - centroide", 
     ylab="k (no. conglomerados)", 
     xlab="h (altura del nodo)", 
     col="grey")
text(espeNEU.CEN$height, 
     nrow(especiesC):2, 
     nrow(especiesC):2, 
     col="red", cex=0.8)

# Número de conglomerados por Ward
plot(espeNEU.WA$height, nrow(especiesC):2, type="S", 
     main="Fusion levels - Cuerda - Ward", 
     ylab="k (no. conglomerados)", 
     xlab="h (altura del nodo)", 
     col="grey")
text(espeNEU.WA$height, 
     nrow(especiesC):2, 
     nrow(especiesC):2, 
     col="red", cex=0.8)
plot(hclust(vegdist(decostand(especiesC, "norm"), "euclidea"), "average"))

# Escogemos un número de conglomerado
k <- 10
library(vegan)
closVC <- cutree(espeNEU.VC, k)
closVL <- cutree(espeNEU.VL, k)
closUP <- cutree(espeNEU.UPGMA, k)
closCE <- cutree(espeNEU.CEN, k)
closWA <- cutree(espeNEU.WA, k)

# Comparar agrupaciones por tablas de contingencia
table(closVC, closVL)
table(closVC, closUP)
table(closVC, closCE)
table(closVC, closWA)

# Obtención de número óptimo de conglomerados por siluetas
asw <- numeric(nrow(especiesC))

# Second, retrieve the asw values and write them into the vector
for (k in 2:(nrow(especiesC)-1)) {
  sil <- silhouette(cutree(espeNEU.WA, k=k), espeNEU)
  asw[k] <- summary(sil)$avg.width
}

kM <- which.max(asw)
dev.new(title="Silhouettes - Ward - k = 2 to n-1")
plot(1:nrow(especiesC), asw, type="h", 
     main="Silhouette-optimal number of clusters, Ward", 
     xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, kM, paste("optimum",kM,sep="\n"), col="red", font=2,
     col.axis="red")
points(kM, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", kM, "\n", 
    "with an average silhouette width of", max(asw), "\n")

### Análisis de perfil de similitud SIMPROF
library(clustsig)

perfil <- simprof(espeNEU, num.expected = 1000,
        method.transform = "squareroot",
        method.distance = "euclidean",
        method.cluster = "average",
        alpha = 0.001)

perfil

dev.new()
simprof.plot(perfil)

perfil <- simprof(especiesC, num.expected = 10,
        method.transform = "squareroot",
        method.distance = "euclidean",
        method.cluster = "ward.D",
        alpha = 0.05)

simprof.plot(perfil)

summary(perfil)

dev.new()
plot(anosim(especiesA, especies$habitat))


rownames(autos) <- paste0(c(1:48),autos$tipo, autos$zona)
autosA <- autos[c(5:38)]

au1 <- decostand(autosA, "norm")
au2 <- vegdist(au1, "euclidea")
au3 <- hclust(au2, "ward.D")
dev.new()
plot(au3)

perfil2 <- simprof(au2, num.expected = 1000,
                  method.transform = "squareroot",
                  method.distance = "euclidean",
                  method.cluster = "ward.D",
                  alpha = 0.05)

simprof.plot(perfil2)

plot(anosim(autosA, autos$zona))


### K medias -------------------------------------
inPr <- read.csv("Datos/promedios.csv", header = T)

especiesA <- especies[c(-1:-3)]
indicesA <- indices[c(5:10)]

especiesC <- rowsum(especiesA, 
                    paste0(especies$punto, especies$habitat))

### Regresar -------------------------------
esN1<-decostand(especiesC, "norm")
esE1 <- vegdist(esN1, "euclidea")
esC1 <- hclust(esE1, "ward.D2")
dev.new()
plot(esC1) ####################################


library(cluster)
k <- 3
esNC <- decostand(especiesC, "norm")
esME <- vegdist(esNC, "euclidea")
esMEC <- hclust(esME, "ward.D2")
mediask <- kmeans(especiesC, centers=3, nstart=100)

esCAS <- cascadeKM(esNC,
                   inf.gr=2, sup.gr=3, 
                   iter=100,
                   criterion="ssi")


esCAS$results
esCAS$partition


esPAM <- pam(especiesC, k=3, diss=TRUE)
summary(esPAM)
esPAMC <- esPAM$clustering
esPAM$silinfo$widths

sil <- silhouette(mediask$cluster, esME)
rownames(sil) <- row.names(especiesC)
dev.new()
par(mfrow=c(1,2))
plot(sil, main="Silueta - k-medias", 
     cex.names=0.8, col=2:(k+1))

plot(silhouette(esPAM), 
     main="Silueta - PAM", cex.names=0.8, 
     col=2:(k+1))


library(factoextra)
dev.new()
fviz_cluster(pam(especiesC, 3),
             ggtheme = theme_classic())

closk <- mediask$partition[,1]
closk <- mediask$cluster

dev.new()
#### Recordar fcuando termine el modo R -----------------------
boxplot(inPr$dba ~mediask$cluster)

hist(inPr$dba)
shapiro.test(inPr$dba)

summary(glm(inPr$dba ~ as.character(mediask$cluster), Gamma))
