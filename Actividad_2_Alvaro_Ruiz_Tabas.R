setwd("/media/alvaro/DISCO_DURO_ART/Transcriptomica_actividad_2/")

### Pregunta 1
# Cargamos diseño experimental (metadatos)
metadata <- read.delim("experimento_GSE167749.tsv", row.names = 1)
par(mar=c(6,6,3,3))
heatmap(table(metadata$disease_state, metadata$GEO))

# Matriz de diseño
grupos <- factor(metadata$disease_state)
diseño <- model.matrix(~0+grupos)
colnames(diseño) <- levels(grupos)
colnames(diseño)[3] <- "infected_treated"
rownames(diseño) <- rownames(metadata)
diseño

### Pregunta 2
# Cargamos datos de recuentos
recuentos <- read.delim("recuento_GSE167749.tsv", row.names = 1)
colnames(recuentos) <- rownames(metadata)
head(recuentos)

# Relacion media-varianza
mean_counts <- apply(recuentos[,6:10], 1, mean)
variance_counts <- apply(recuentos[,6:10], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  ggtitle("Relacion media-varianza en muestras \n enfermas no tratadas") +
  xlab("Media") + ylab("Varianza")

# Genes en la matriz de recuentos
dim(recuentos)
length(grep("ENSMUS", rownames(recuentos), fixed=TRUE)) # Genes de ratón
length(grep("ENSMUS", rownames(recuentos), fixed=TRUE, invert = TRUE)) # Genes de COVID

# Eliminacion de genes de expreión nula, 4 y 10.
library(statmod)
library(edgeR)
recuentos_filtrado_0 <- recuentos[rowSums(recuentos) > 0,]
nrow(recuentos_filtrado_0)
par(mar=c(6,6,3,3))
hist(aveLogCPM(recuentos_filtrado_0), 
     main= 'Dataset tras eliminar genes \n de expresión nula',
     xlab= "Media logCPM",
     ylab="Número de genes",
     xlim= c(-5,15),
     ylim= c(0,10000),
     col = "blue"
)

recuentos_filtrado_4 <- recuentos[rowSums(recuentos) > 4,]
nrow(recuentos_filtrado_4)
par(mar=c(6,6,3,3))
hist(aveLogCPM(recuentos_filtrado_4), 
     main= 'Dataset tras eliminar genes \n de expresión > 4',
     xlab= "Media logCPM",
     ylab="Número de genes",
     xlim= c(-5,15),
     ylim= c(0,10000),
     col = "yellow"
)

recuentos_filtrado_10 <- recuentos[rowSums(recuentos) > 10,]
nrow(recuentos_filtrado_10)
par(mar=c(6,6,3,3))
hist(aveLogCPM(recuentos_filtrado_10), 
     main= 'Dataset tras eliminar genes \n de expresión > 10',
     xlab= "Media logCPM",
     ylab="Número de genes",
     xlim= c(-5,15),
     ylim= c(0,10000),
     col = "red"
)

### Pregunta 3
# Objeto DGElist y normalizacion
y <- DGEList(recuentos_filtrado_4, group = grupos)
y_norm <- calcNormFactors(y, method = "TMM")

# Factores de normalización por muestra
y_norm$samples

# Boxplot de los contajes sin y con normalizacion

cpms_NO_norm <- cpm(y, normalized.lib.sizes = F, log=T)
boxplot(cpms_NO_norm, xlab="muestras", ylab="Log2 CPM",las=2, col = "darkolivegreen1")
abline(h=median(cpms_NO_norm),col="red")
title("Librerías NO normalizadas(logCPMs)")

cpms_norm <- cpm(y_norm, normalized.lib.sizes = T, log=T)
boxplot(cpms_norm, xlab="muestras", ylab="Log2 CPM",las=2, col = "aquamarine")
abline(h=median(cpms_norm),col="red")
title("Librerías normalizadas(logCPMs)")

# Grafico MDS
plotMDS(y_norm, ylab = "Dimensión 2", xlab = "Dimensión 1", 
        col = as.numeric(y_norm$samples$group), 
        main = "Gráfico MDS")

### Pregunta 4
# Calculo de las diapersiones
y_norm <- estimateDisp(y_norm, diseño, robust = TRUE)

# Trended y tagwise dispersion para genes concretos
y_norm$trended.dispersion[which(rownames(y_norm$counts) == "ENSMUSG00000078354")]
y_norm$tagwise.dispersion[which(rownames(y_norm$counts) == "ENSMUSG00000078354")]

y_norm$trended.dispersion[which(rownames(y_norm$counts) == "ENSMUSG00000049176")]    
y_norm$tagwise.dispersion[which(rownames(y_norm$counts) == "ENSMUSG00000049176")]

y_norm$trended.dispersion[which(rownames(y_norm$counts) == "55060280")]
y_norm$tagwise.dispersion[which(rownames(y_norm$counts) == "55060280")]

# Ajuste GLM
fit <- glmQLFit(y_norm, diseño, robust = TRUE)

# Coeficientes para los grupos experimentales en genes concretos
fit$coefficients["ENSMUSG00000078354",]
fit$coefficients["ENSMUSG00000049176",]
fit$coefficients["55060280",]


### Pregunta 5
# Test Infectado vs Sano + corrección con Toptags
Infected_vs_Healthy_contrast <- makeContrasts(infected-healthy, levels = diseño)
Infected_vs_Healthy <- glmQLFTest(fit, contrast = Infected_vs_Healthy_contrast)
Infected_vs_Healthy_DGE <- topTags(Infected_vs_Healthy, n=Inf)

# Test Tratado vs Infectado + corrección con TopTags
Treated_vs_Infected_contrast <- makeContrasts(infected_treated-infected, levels = diseño)
Treated_vs_Infected <- glmQLFTest(fit, contrast = Treated_vs_Infected_contrast)
Treated_vs_Infected_DGE <- topTags(Treated_vs_Infected, n=Inf)

# Histogramas p-value y FDR Infected vs Healthy
par(mar=c(6,6,3,3))
hist(Infected_vs_Healthy_DGE$table$PValue, 
     main= 'Valores de p-value Infectado vs Sano',
     xlab= "Valores de p-value",
     ylab="Frecuencia",
     xlim= c(0.001, 1), 
     ylim= c(0,10000), 
     col = "red", 
     breaks = 40)

par(mar=c(6,6,3,3))
hist(Infected_vs_Healthy_DGE$table$FDR, 
     main= 'Valores de FDR Infectado vs Sano',
     xlab= "Valores de FDR",
     ylab="Frecuencia",
     xlim= c(0.001, 1), 
     ylim= c(0,10000), 
     col = "blue", 
     breaks = 40)

# Histograma p-value y FDR Tratado vs Infectado
par(mar=c(6,6,3,3))
hist(Treated_vs_Infected_DGE$table$PValue, 
     main= 'Valores de p-value Tratado vs Infectado',
     xlab= "Valores de p-value",
     ylab="Frecuencia",
     xlim= c(0.001, 1), 
     ylim= c(0, 10000), 
     col = "red", 
     breaks = 40)

par(mar=c(6,6,3,3))
hist(Treated_vs_Infected_DGE$table$FDR, 
     main= 'Valores de FDR Tratado vs Infectado',
     xlab= "Valores de FDR",
     ylab="Frecuencia",
     xlim= c(0.001, 1), 
     ylim= c(0,10000), 
     col = "blue", 
     breaks = 40)

# DEG con FDR < 0.05 en Infected vs Healthy + plotMD + FDR < 0.01 
nrow(Infected_vs_Healthy_DGE$table)
Infected_vs_Healthy_DGE_0.05 <- Infected_vs_Healthy_DGE$table[Infected_vs_Healthy_DGE$table$FDR <= 0.05,]
nrow(Infected_vs_Healthy_DGE_0.05) # DEG en Infected vs healty con FDR <= 0.05

is.de_Infected_vs_Healthy <- decideTests(Infected_vs_Healthy)
# sum(is.de_Infected_vs_Healthy == 1/-1) > Ver DEG considerados por decideTest. Umbral de FDR <= 0.05
plotMD(Infected_vs_Healthy, status = is.de_Infected_vs_Healthy, cex=0.5, main = "MDplot Infectado vs Sano")
abline(h = c(-1,1), col = "black")

Infected_vs_Healthy_DGE_0.01 <- Infected_vs_Healthy_DGE$table[Infected_vs_Healthy_DGE$table$FDR <= 0.01,]
nrow(Infected_vs_Healthy_DGE_0.01)

# DEG con FDR < 0.05 en Treated vs Infected + plotMD + FDR < 0.01
nrow(Treated_vs_Infected_DGE$table)
Treated_vs_Infected_DGE_0.05 <- Treated_vs_Infected_DGE$table[Treated_vs_Infected_DGE$table$FDR <= 0.05,]
nrow(Treated_vs_Infected_DGE_0.05) # DEG en Infected vs healty con FDR <= 0.05

is.de_Treated_vs_Infected <- decideTests(Treated_vs_Infected)
# sum(is.de_Treated_vs_Infected == 1/-1) > Ver DEG considerados por decideTest. Umbral de FDR <= 0.05
plotMD(Treated_vs_Infected, status = is.de_Treated_vs_Infected, cex=0.5, main = "MDplot Tratado vs Infectado", 
       legend = "bottomright")
abline(h = c(-1,1), col = "black")

Treated_vs_Infected_DGE_0.01 <- Treated_vs_Infected_DGE$table[Treated_vs_Infected_DGE$table$FDR <= 0.01,]
nrow(Treated_vs_Infected_DGE_0.01)

# Volcano plot Infectado vs Sano con filtro de FDR y logFC
library(devtools)
library(EnhancedVolcano)
library(tidyverse)

Infected_vs_Healthy_DGE_0.05_logFC_1 <- Infected_vs_Healthy_DGE_0.05[abs(Infected_vs_Healthy_DGE_0.05$logFC)>=1,]

nrow(Infected_vs_Healthy_DGE_0.05_logFC_1) # Genes que pasan ambos filtros
sum(Infected_vs_Healthy_DGE_0.05_logFC_1$logFC > 0) # Genes sobreexpresados
sum(Infected_vs_Healthy_DGE_0.05_logFC_1$logFC < 0) # Genes infraexpresados


Infected_vs_Healthy_DGE_0.05_logFC_1 <- Infected_vs_Healthy_DGE_0.05_logFC_1 %>%
  arrange(desc(abs(Infected_vs_Healthy_DGE_0.05_logFC_1$logFC)))
top3_logFC <- rownames(Infected_vs_Healthy_DGE_0.05_logFC_1[1:3,]) # Genes con mayor LogFC 

Infected_vs_Healthy_DGE_0.05_logFC_1 <- Infected_vs_Healthy_DGE_0.05_logFC_1 %>%
  arrange(Infected_vs_Healthy_DGE_0.05_logFC_1$FDR)
top3_FDR <- rownames(Infected_vs_Healthy_DGE_0.05_logFC_1[1:3,]) # Genes con menor FDR

top_genes <- c(top3_FDR, top3_logFC) # Top genes segun FDR y log FC

EnhancedVolcano(Infected_vs_Healthy_DGE$table, lab = rownames(Infected_vs_Healthy_DGE$table),
                selectLab = top_genes,
                boxedLabels = TRUE, drawConnectors = TRUE, subtitle = "",
                x = 'logFC', y = 'FDR', pCutoff = 0.05, FCcutoff = 1, title = "Volcano plot Infectado vs Sano")

# Volcano plot Tratado vs infectado con filtro de FDR y logFC
Treated_vs_Infected_DGE_0.05_logFC_1 <- Treated_vs_Infected_DGE_0.05[abs(Treated_vs_Infected_DGE_0.05$logFC)>=1,]

nrow(Treated_vs_Infected_DGE_0.05_logFC_1) # Genes que pasan ambos filtros
sum(Treated_vs_Infected_DGE_0.05_logFC_1$logFC > 0) # Genes sobreexpresados
sum(Treated_vs_Infected_DGE_0.05_logFC_1$logFC < 0) # Genes infraexpresados


Treated_vs_Infected_DGE_0.05_logFC_1 <- Treated_vs_Infected_DGE_0.05_logFC_1 %>%
  arrange(desc(abs(Treated_vs_Infected_DGE_0.05_logFC_1$logFC)))
top3_logFC <- rownames(Treated_vs_Infected_DGE_0.05_logFC_1[1:3,]) # Genes con mayor LogFC 

Treated_vs_Infected_DGE_0.05_logFC_1 <- Treated_vs_Infected_DGE_0.05_logFC_1 %>%
  arrange(Treated_vs_Infected_DGE_0.05_logFC_1$FDR)
top3_FDR <- rownames(Treated_vs_Infected_DGE_0.05_logFC_1[1:3,]) # Genes con menor FDR

top_genes <- c(top3_FDR, top3_logFC) # Top genes segun FDR y log FC


EnhancedVolcano(Treated_vs_Infected_DGE$table, lab = rownames(Treated_vs_Infected_DGE$table),
                selectLab = top_genes,
                boxedLabels = TRUE, drawConnectors = TRUE, subtitle = "",
                x = 'logFC', y = 'FDR', pCutoff = 0.05, FCcutoff = 1, title = "Volcano plot Tratado vs Infectado")

# Diagrama de Venn y Upset plot entre grupos comparados
Infected_vs_Healthy_UP <- rownames(Infected_vs_Healthy_DGE_0.05_logFC_1
                                   [Infected_vs_Healthy_DGE_0.05_logFC_1$logFC > 0,])
Infected_vs_Healthy_DOWN <- rownames(Infected_vs_Healthy_DGE_0.05_logFC_1
                                     [Infected_vs_Healthy_DGE_0.05_logFC_1$logFC < 0,])
Treated_vs_Infected_UP <- rownames(Treated_vs_Infected_DGE_0.05_logFC_1
                                   [Treated_vs_Infected_DGE_0.05_logFC_1$logFC > 0,])
Treated_vs_Infected_DOWN <- rownames(Treated_vs_Infected_DGE_0.05_logFC_1
                                     [Treated_vs_Infected_DGE_0.05_logFC_1$logFC < 0,])

library(venn)
venn(list("I_vs_H_up"=Infected_vs_Healthy_UP, "I_vs_H_down"=Infected_vs_Healthy_DOWN, 
          "T_vs_I_up"=Treated_vs_Infected_UP, "T_vs_I_down"=Treated_vs_Infected_DOWN ),
     zcolor = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

library(UpSetR)
InputList <- list("I_vs_H_up"=Infected_vs_Healthy_UP, "I_vs_H_down"=Infected_vs_Healthy_DOWN, 
                  "T_vs_I_up"=Treated_vs_Infected_UP, "T_vs_I_down"=Treated_vs_Infected_DOWN )

upset(fromList(InputList), order.by = "freq", number.angles = 30, 
      point.size = 3.5, line.size = 2, mainbar.y.label = "Groups Intersections",
      sets.x.label = "Genes per group")


# Pregunta 6
# Transformamos a z-scores a partir d elos counts normalizados 
library(matrixStats)
cpms <- cpm(y_norm$counts, log = TRUE)
colnames(cpms) <- rownames(metadata)
all_DGE_filt <- intersect(row.names(Infected_vs_Healthy_DGE_0.05_logFC_1),row.names(Treated_vs_Infected_DGE_0.05_logFC_1))
matrix_DGE_cpm <- cpms[all_DGE_filt,]
zscore_fun <- function(x){(x - mean(x))/sd(x)}
matrix_DGE_zscore <- t(apply(matrix_DGE_cpm, 1, zscore_fun))
head(rowSds(matrix_DGE_zscore))

par(mar=c(6,6,3,3))
d <- density(matrix_DGE_cpm)
plot(d, main = "Gráfico de densidad CPMs") # Density plot CPMs
polygon(d, col="red", border="blue")

d2 <- density(matrix_DGE_zscore)
plot(d2, main = "Gráfico de densidad Zscore") # Density plot Zscores
polygon(d2, col="blue", border="red")

# Heatmap
library(pheatmap)
# Calculo distancia euclidianas
heatmap <- pheatmap(
  matrix_DGE_cpm, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap clusterizado por metodo Euclidean",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "deepskyblue",
    "black",
    "yellow"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)

# Calculo distancias Manhattan
heatmap <- pheatmap(
  matrix_DGE_cpm, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "manhattan",
  clustering_distance_cols = "manhattan",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap clusterizado por metodo Manhattan",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "red",
    "white",
    "darkgreen"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)

# Calculo distancias maximum
heatmap <- pheatmap(
  matrix_DGE_cpm, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "maximum",
  clustering_distance_cols = "maximum",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap clusterizado por metodo Maximum",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "darkblue",
    "white",
    "red4"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)

# Heatmap solo con genes de covid

matrix_DGE_cpm_COVID <- matrix_DGE_cpm[grep("ENSMUS", rownames(matrix_DGE_cpm), 
                                            fixed=TRUE, invert = TRUE),]
heatmap <- pheatmap(
  matrix_DGE_cpm_COVID, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap genes de COVID con linkage complete",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "deepskyblue",
    "black",
    "yellow"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)

heatmap <- pheatmap(
  matrix_DGE_cpm_COVID, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap genes de COVID con linkage average",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "red",
    "white",
    "darkgreen"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)


heatmap <- pheatmap(
  matrix_DGE_cpm_COVID, # Metemos los datos nosmalizados SOLO, no escalados (se lo indicamos abajo)
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "single",
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Heatmap genes de COVID con linkage single",
  treeheight_row = 0,
  #treeheight_col = 0,
  colorRampPalette(c(
    "darkblue",
    "white",
    "red4"))(25),
  scale = "row" # Scale values in the direction of genes (rows)
  # Si cambiamos por none se ve sin escalar
)


# Pregunta 7
require(org.Mm.eg.db)

base <- as.data.frame(org.Mm.egENSEMBL2EG)

funer <- base[base$ensembl_id %in% Infected_vs_Healthy_UP, ]
go <- goana(funer$gene_id, species="Mm", geneid= "ENTREZID") # Busqueda de go y estudio de su frecuencia 
topGO(go, ontology="BP", number=3)
topGO(go, ontology="MF", number=3)

funer <- base[base$ensembl_id %in% Infected_vs_Healthy_DOWN, ]
go <- goana(funer$gene_id, species="Mm", geneid= "ENTREZID") # Busqueda de go y estudio de su frecuencia 
topGO(go, ontology="BP", number=3)
topGO(go, ontology="MF", number=3)

funer <- base[base$ensembl_id %in% Treated_vs_Infected_UP, ]
go <- goana(funer$gene_id, species="Mm", geneid= "ENTREZID") # Busqueda de go y estudio de su frecuencia 
topGO(go, ontology="BP", number=3)
topGO(go, ontology="MF", number=3)

funer <- base[base$ensembl_id %in% Treated_vs_Infected_DOWN, ]
go <- goana(funer$gene_id, species="Mm", geneid= "ENTREZID") # Busqueda de go y estudio de su frecuencia 
topGO(go, ontology="BP", number=3)
topGO(go, ontology="MF", number=3)



write.csv(Infected_vs_Healthy_UP, "Infected_vs_Healthy_UP.csv", row.names = FALSE, col.names = FALSE)
