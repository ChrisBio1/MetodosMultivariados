especiesC <- rowsum(especiesA, paste0(especies$punto, especies$habitat))

### Análisis de conglomerados --------------------------------------
# Método de aglomeramiento del vecino más cercano
espN <- decostand(especiesC, "normalize")
espeNEU <- vegdist(espN, "euc")
espeNEU.VC <- hclust(espeNEU, method="single")

plot(espeNEU.VC)
plot(hclust(vegdist(decostand(especiesC, "normalize"),"euclidean"), "single"))

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
m1 <- head(autos, 8)[c(5:38)]

plot(hclust(vegdist(decostand(m1, "norm") ,"euclidean"), "ward.D"))
###

### Correlación cofenética -------------------------------------
# Vecino más cercano
espeCofC <- cophenetic(espeNEU.VC)
cofC <- cor(espeNEU, espeCofC)

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

# Número de conglomerados por vecino más cercano
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

# Escogemos un número de conglomerado
k <- 5

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
