especies <- read.csv("Datos/especies.csv", header = T)
indices <- read.csv("Datos/indices.csv", header = T)
ubicacion <- read.csv("Datos/ubicacion.csv", header = T)

capitaliza <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


### Explorar Bases de datos -----------------------------------------------

### Etiquetas de las columnas de las bases de datos
names(especies)
names(indices)

### Dimensiones de las bases de datos
dim(indices)
dim(especies)

### Observar las primeras 6 filas de cada columna
head(indices)

### Observar las últimas 6 filas de cada columna
tail(indices)

### Extraer filas o columnas de bases de datos
### Mediante []
# Primera columna
indices[1]
indices[-1]
# Más de una columna
indices[c(1:3)]
indices[c(-1:-3)]
indices[c(-1:-3,-6)]

# Primera fila
indices[1,]
indices[-1,]
# Más de dos filas
indices[c(1:3),]
indices[c(-1:-3),]
indices[c(-1:-3, -7:-9),]
# Por grupo
especies[especies$habitat == "cultivo",]


### Mediante $
# La variable hábitat
indices$habitat


### Unir bases de datos
# Por columnas
cbind(indices, especies)


# A partir de un criterio
merge(indices, especies, by = "habitat")

# Crear criterio
id <- paste0(indices$dia, indices$habitat, indices$punto)

# Unir criterio a las dos bases de datos y renombrar
especies2 <- cbind(id, especiesA)
indices2 <- cbind(id, indices)

# Unir bases mediante criterio agregado
merge(indices2, especies2, by = "id")

# Transponer matriz
t(especies)

### Análisis exploratorio -------------------------------------
# Medidas de tendencia central
tapply(indices$dba, indices$habitat, mean)
tapply(indices$dba, indices$habitat, median)

# Pruebas de normalidad para todas las variables de índices
lapply(indices[c(-1:-4)], function(x){
  shapiro.test(x)
})

# Histogramas de algunas especies
par(mfrow = c(1,3))
lapply(especies[c(30:32)], function(x){
  hist(x, xlab = "", main = names(x))
})

### Descriptivos de las especies 
especiesA <- especies[c(-1:-3)]

### Valores mínimo y máximo de abundancias por punto de muestreo
range(especiesA)

### Frecuencia de abundancias
ab <- table(unlist(especiesA))

### Gráfica de barras de frecuencias de abundancias
barplot(ab, las=1, xlab="Abundancias", ylab="Frecuencia", col=gray(5:0/5))

### Número de ausencias
sum(especiesA==0)

### Proporción de ausencias
sum(especiesA==0)/(nrow(especiesA)*ncol(especiesA))

# Riqueza por hábitat
s <- rowSums(ifelse(rowsum(especiesA, 
                           especies$habitat) > 0, 1, 0))

# Ausencias por hábitat
si <- rowSums(ifelse(rowsum(especiesA, 
                      especies$habitat) == 0, 1, 0))

si/(s+si)*100

### 
### Presencias totales
ab <- rowSums(rowsum(ifelse(especiesA > 0,1, 0), 
                     especies$habitat))

### Ausencias totales
au <- rowSums(rowsum(ifelse(especiesA == 0,1, 0), 
                     especies$habitat))

### Porcentaje de ausencias
mean(au/(ab+au)*100)

# Riqueza por punto de muestreo
sm <- rowSums(especiesA > 0)
data.frame(especies$habitat, sm)

# Abundancia de especies por hábitat
abunni <- rowsum(especiesA, especies$habitat)
# Número total de individuos en cada hábitat
abunN <- colSums(t(abunni))

# Cálculo de pi
Api <- abunni/abunN

# Alfa de Simpson
rowSums(Api^2)

# Alfa de Shannon-Wiener (H)
H1 <- rowSums(-(Api*log(Api)), na.rm = T)

# Equidad de Pielou por hábitat
H1/log(s)

N0 <- rowSums(abunni > 0)                # Species richness
H <- vegan::diversity(abunni)            # Shannon entropy
N1 <- exp(H)                   # Shannon diversity (number of abundant species)
N2 <- vegan::diversity(abunni, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div <- data.frame(N0, H, N1, N2, E10, E20, J))


head(indices)

### Rango abundancia 
# Para un grupo
a1 <- rowSums(t(especies[especies$habitat == "cultivo",][c(-1:-3)]))

par(bty = "n")
par(family = "serif")
plot(c(1:length(a1)), sort(a1, decreasing = T), 
     type = "l",
     xlim = c(-8, length(a1)),
     ylim = c(0, max(a1)+40),
     ylab = "",
     xlab = "",
     las =1,
     xaxt = "n")
points(c(1:length(a1)), sort(a1, decreasing = T),
       pch = 16)
n1 <- head(a1,6)
text(head(c(1:length(a1)), 6),
     head(sort(a1, decreasing = T), 6), 
     gsub("[.]", " ", names(n1)),
     pos = c(4,2), 
     cex = 1,
     srt = 45,
     font = 3)
text(40, 50, 
     c(paste0("S = ",length(a1))),
     cex = 2)
mtext(capitaliza(levels(especies$habitat)[1]),
      1, 1, cex = 1.5)
mtext("Abundancia", 2, 2.5, cex = 1.5)


### Para más de un grupo
rangos <- list()
for(i in levels(especies$habitat)){
  abundancia <- t(sort(rowsum(especies[c(-1,-2, -3)], 
                              especies$habitat)[i,],
                       decreasing = T))
  abundancia <- na.omit(ifelse(abundancia == 0, NA, abundancia))
  ejeX <- c(1:length(abundancia))
  especie <- gsub("[.]", " ", row.names(abundancia))
  rangos[[i]] <- data.frame(especie, abundancia, ejeX)
}

rangos

par(mfrow = c(1,3))
par(family = "serif")
par(bty = "n")
lapply(rangos, function(x){
  plot(x[[3]], x[[2]],
       las = 1, 
       xlab = "", ylab = "",
       ylim = c(0, 125),
       xlim = c(-30, 50),
       col = "cornflowerblue",
       type= "l",
       xaxt = "n",
       lwd = 2, 
       cex.axis = 1.5)
  points(x[[3]], x[[2]], 
         pch = 16,
         cex = 1.5)
  x2 <- head(x,6)
  text(x2[[3]],
       x2[[2]], 
       x2[[1]],
       pos = c(4,2), 
       cex = 2.7,
       srt = 45,
       font = 3)
  text(40, 50, 
       c(paste0("S = ",length(x[[3]]))),
       cex = 2.5)
})
mtext("Abundancia", 2, 83, cex = 1.5)
# mtext("Abundancia", 2, 40, cex = 1.5)
mtext(capitaliza(levels(especies$habitat)), 
      rep(1, 3), 2, 0, c(-190, -80, 15), 
      cex = 2)

### Datos espaciales
ub2 <- ubicacion$habitat
levels(ub2) <- c(rgb(1,0,0, alpha = 0.5),
                 rgb(0,1,0, alpha = 0.5), 
                 rgb(0,0,1, alpha = 0.5))
ubicacion$habitat
plot(ubicacion$longitud, ubicacion$latitud, 
     asp = 1, 
     pch = 16, 
     cex = log(sm),
     col = as.character(ub2),
     xlab = "",
     ylab = "")
text(ubicacion$longitud, 
      ubicacion$latitud,
     ubicacion$punto)
text(ubicacion$longitud[c(25, 33, 49)], 
      ubicacion$latitud[c(25, 33, 49)],
     ubicacion$habitat[c(25, 33, 49)],
     cex = 2)
mtext("Longitud", 1, 2, cex = 1.5)
mtext("Latitud", 2, 2, cex = 1.5)

### Uso de "Vegan" ---------------------------------------------------
library(vegan)
### Transformación y estandarización de datos de especies ------------
decostand(especiesA, "pa")

### Escala de abundancias por su división por el valor máximo de cada especie
escala <- decostand(especiesA, "max")

# Obtener el valor máximo de cada columna (especie)
apply(escala, 2, max)

### Escalar abundancias dividiendo por el total de especies
totalE <- decostand(especiesA, "total", MARGIN=2)
totalE[1:5,2:4]

colSums(totalE)


### Transformaciones útiles para PCA o RDA
### Transformación de cuerda
cuerda <- decostand(especiesA, "normalize") # default MARGIN=1
cuerda[1:5,2:4]

colSums(cuerda)


### Verificar normalización de vectores de filas
norm <- function(x) sqrt(x%*%x)
apply(cuerda, 1, norm)
apply(totalE, 1, norm)


### Raíz cuadrada de las abundancias relativas por sitio
hel <- decostand(especiesA, "hellinger")
hel[1:5,2:4]

apply(hel, 1, norm)

### Estandarización por filas y columnas
### Tranformación Chi cuadrada
chi <- decostand(especiesA, "chi.square")
chi[1:5,2:4]

apply(chi, 1, norm)


### Estandarización de Wisconsin
wis <- wisconsin(especiesA)
wis[1:5,2:4]

apply(wis, 1, norm)


### Diferenciast
par(mfrow= c(2,2))
par(bty = "n")
boxplot(cbind(especies[4], sqrt(especies[4]), log1p(especies[4])),
        las=1, main="Transformación simple",
        names=c("Frecuencia", "sqrt", "log"), 
        col="lavender", 
        horizontal = F, 
        border = rgb(0, 0, 1, alpha = 0.5), 
        boxwex = .5)

boxplot(cbind(escala[1], totalE[1]),
        las=1, main="Estandarización por especies",
        names=c("Máximo", "Total"),
        col="lavender", 
        horizontal = F, 
        border = rgb(0, 0, 1, alpha = 0.5), 
        boxwex = .5)

boxplot(cbind(cuerda[1], hel[1], totalE[1]),
        las=1, main="Estandarización por sitios",
        names=c("Cuerda", "Hellinger", "Total"),        
        col="lavender", 
        horizontal = F, 
        border = rgb(0, 0, 1, alpha = 0.5), 
        boxwex = .5)

boxplot(cbind(chi[1], wis[1]),
        las=1, main="Doble estandarización",
        names=c("X²", "Wisconsin"),
        col="lavender", 
        horizontal = F, 
        border = rgb(0, 0, 1, alpha = 0.5), 
        boxwex = .5)

### MODO Q -----------------------------------------------------------
# Medidas de similitud y distancias para datos cuantitavios
esB <- vegdist(especiesA)	# method="bray" (default)
head(esB)

# Percentage difference (Bray-Curtis) dissimilarity matrix
# on log-transformed abundances
esL <- vegdist(log1p(especiesA))
head(esL)

# Chord distance matrix
esN <- decostand(especiesA, "nor")
esC <- dist(esN)
head(esC)

# Hellinger distance matrix
esH <- decostand(especiesA, "hel")
esH2 <- dist(esH)
head(esH2)

#### Modo Q medidas de disimilitud para datos binarios ####
# *********************************************

# Disimilitud de Jaccard con vegdist
especiesB <- ifelse(especiesA > 0, 1, 0)

esJ <- vegdist(especiesB, "jac", binary=TRUE)
head(esJ)
head(sqrt(esJ))

# Disimilitud de Jaccard con dist()
esB <- dist(especiesB, "binary")
head(esB)

# Disimilitud de Sorensen con vegdist()
esS <- vegdist(especiesB, binary=TRUE)
head(esS)
head(sqrt(esS))

### ------
source("Scripts/coldiss.R")
source("Scripts/panelutils.R")

# Percentage difference (Bray-Curtis) dissimilarity matrix on raw species abundance data
dev.new(title="Bray Curtis", width=10, height=5)
coldiss(esB, byrank=FALSE, diag=TRUE)

# Matriz de distancias por logarítmo de Bray Curtis
dev.new(title="Bray-Curtis, ln(y+1)", width=10, height=5)
coldiss(esL, byrank=FALSE, diag=TRUE)

# Matriz de distancia de cuerda
dev.new(title="Cuerda", width=10, height=5)
coldiss(esC, byrank=FALSE, diag=TRUE)

# Matriz de distancias de Hellinger
# dev.new(title="Hellinger", width=10, height=5)
# coldiss(esH, byrank=FALSE, diag=TRUE)

# Matriz de distancias de Jaccard
dev.new(title="Jaccard", width=10, height=5)
coldiss(esJ, byrank=FALSE, diag=TRUE)


