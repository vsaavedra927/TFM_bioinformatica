rm(list=ls())


# CARGAMOS LIBRERIAS

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

# Semilla de aleatorizacion
set.seed(1989)


# Para informacion detallada de DESeq2, ejecutar esta linea:
#vignette("DESeq2")


#### ---- CARGA DE DATOS ---- ####

# NOTAS: 
# DATOS: Es importante usar CONTEOS NO NORMALIZADOS como input para DESeq2 (no normalizados, no transformados)
# DATOS vs METADATOS: los nombres de muestras deben encontrarse como nombres de columnas en los DATOS y como nombres de filas en los METADATOS, y en el mismo orden.

# CARGAMOS RUTA Y NOMBRES DE ARCHIVO

# Indicar codigo de experimento y descripcion (p.ej. con autor). Para rutas de archivo
gse_nr <- "GSE87426"
gse_description <- "GSE87426_Ferguson2016"

# Marcamos la ruta y guardamos nombres de archivo en una lista
path_ubuntu <- paste0("//wsl.localhost/Ubuntu/home/victor/UNIR/TFM/Datasets/", gse_description, "/")
setwd(path_ubuntu)
file_names <- list.files(path = path_ubuntu, pattern = "\\.tsv$")
file_names
data_file <- file_names[1]

# Marcamos la ruta para guardar resultados. Crearla si no existe
path_results <- paste0("//wsl.localhost/Ubuntu/home/victor/UNIR/TFM/Results/Results_", gse_description, "/")
if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}


# CREAMOS DATAFRAME DE EXPRESION

# Existe un solo archivo para abrir en un dataframe

# Construimos ruta del archivo
file_path <- paste0(path_ubuntu, data_file)

# Leemos el archivo con los datos
exp_df_all <- as.data.frame(read_tsv(file_path))

# Exploramos rapidamente las primeras filas y la estructura del dataframe de expresion
head(exp_df_all)
str(exp_df_all)
# Cambiaremos los codigos de genes usando las anotaciones en el archivo de anotaciones


# Abrimos archivo de anotaciones
annotations_file <- file_names[2]
annot_path <- paste0(path_ubuntu, annotations_file)
annotations_df <- read_tsv(annot_path)

# Comprobar que genes en df de datos y de anotaciones coinciden
print(paste0("Coinciden genes: ", all(annotations_df$GeneID == exp_df_all$GeneID)))
annotations_df$GeneID == exp_df_all$GeneID

# Ponemos nombres de genes como nombres de fila en df de datos y eliminamos columna de codigos de genes
exp_df_all$GeneID <- NULL
rownames(exp_df_all) <- make.unique(annotations_df$Symbol)



# CREAMOS DATAFRAME DE METADATOS

# Acceso a metadatos usando getGEO
gse <- getGEO(gse_nr, GSEMatrix = TRUE)
names(gse)
names(pData(gse[[1]]))
metadata_df_all <- pData(phenoData(gse[[1]]))

# Adaptamos la columna title (el numero de paciente) a una columna subject
metadata_df_all$subject <- sub("^(Subject_\\d+).*", "\\1", metadata_df_all$title)

# Simplificamos las columna treatment y time (nombre, eliminar espacios en blanco, aclarar nombres)
metadata_df_all <- metadata_df_all %>%
  rename(treatment = `treatment:ch1`) %>%
  mutate(treatment = ifelse(treatment == "n-3 PUFA", "Omega3", treatment))

metadata_df_all <- metadata_df_all %>%
  rename(time = `time:ch1`) %>%
  mutate(
    time = case_when(
      time == "Pre-treatment" ~ "T0_Basal",
      time == "Post-treatment Pre-LPS" ~ "T1_Suplementado",
      time == "Post-LPS" ~ "T2_Posttratamiento",
      TRUE ~ time   # Deja otros valores como estan
    )
  )

# CONVERTIR COLUMNAS DE INTERES A FACTOR

# Convertimos variables de interes a factor y marcamos nivel de referencia
metadata_df_all$treatment <- factor(metadata_df_all$treatment) 
metadata_df_all$treatment <- relevel(metadata_df_all$treatment, ref = "Placebo")
levels(metadata_df_all$treatment)

metadata_df_all$time <- factor(metadata_df_all$time) 
levels(metadata_df_all$time)

metadata_df_all$subject <- factor(metadata_df_all$subject) 
levels(metadata_df_all$subject)



# AJUSTE DE DATOS Y METADATOS

# Nos aseguramos de que los nombres de muestras estan en el mismo orden como columnas de datos y filas de metadatos
print(paste0("Coinciden los nombres de todas las muestras: ", all(row.names(metadata_df_all) == colnames(exp_df_all))))
row.names(metadata_df_all) == colnames(exp_df_all)


# DIVISION DEL DATASET EN SUB-DATASETS

# Numeros de muestra para cada timepoint
samples_t0 <- paste0("GSM", seq(2331195, 2331208))
samples_t1 <- paste0("GSM", seq(2331209, 2331222))
samples_t2 <- paste0("GSM", seq(2331223, 2331236))


# Metadata de cada tejido
metadata_df_t0 <- metadata_df_all[samples_t0, ]
metadata_df_t1 <- metadata_df_all[samples_t1, ]
metadata_df_t2 <- metadata_df_all[samples_t2, ]

# Matriz de conteo de cada tejido
exp_df_t0 <- exp_df_all[,samples_t0]
exp_df_t1 <- exp_df_all[,samples_t1]
exp_df_t2 <- exp_df_all[,samples_t2]




### --- OBJETO DESEQDATASET

# OBJETOS DESEQDATASET
# DESeqDataSet de los subsets (distintos timepoints)
dds_t0 <- DESeqDataSetFromMatrix(countData =  exp_df_t0,
                                  colData =  metadata_df_t0,
                                  design = ~ treatment)
dds_t1 <- DESeqDataSetFromMatrix(countData =  exp_df_t1,
                                    colData =  metadata_df_t1,
                                    design = ~ treatment)
dds_t2 <- DESeqDataSetFromMatrix(countData =  exp_df_t2,
                                    colData =  metadata_df_t2,
                                    design = ~ treatment)
# Distintos DESeqDataSet del dataset global. 
# Recordar que subject y treatment en Ferguson2016 son colineales.
dds_all_timebysubject <- DESeqDataSetFromMatrix(countData =  exp_df_all,
                                    colData =  metadata_df_all,
                                    design = ~ subject + time)
dds_all_timetreatment<- DESeqDataSetFromMatrix(countData =  exp_df_all,
                                  colData =  metadata_df_all,
                                  design = ~ time + treatment + time:treatment)


# NORMALIZACION DE LOS CONTEOS

# Determinar los size factors a usar en la normalizacion
dds_t0 <- estimateSizeFactors(dds_t0)
dds_t1 <- estimateSizeFactors(dds_t1)
dds_t2 <- estimateSizeFactors(dds_t2)
dds_all_timebysubject <- estimateSizeFactors(dds_all_timebysubject)
dds_all_timetreatment <- estimateSizeFactors(dds_all_timetreatment)

# Extraer el dataframe de conteos normalizados
exp_df_t0_normalizado <- counts(dds_t0, normalized=TRUE)
exp_df_t1_normalizado <- counts(dds_t1, normalized=TRUE)
exp_df_t2_normalizado <- counts(dds_t2, normalized=TRUE)
exp_df_all_normalizado_timebysubject <- counts(dds_all_timebysubject, normalized=TRUE)
exp_df_all_normalizado_timetreatment <- counts(dds_all_timetreatment, normalized=TRUE)



#### ---- EXPLORACION DE LOS DATOS ---- ####

# CLUSTERIZACION NO SUPERVISADA (heatmap y PCA)
# 1. Hacer un analisis de calidad de las muestras (se agrupan muestras similares entre si como esperamos? hay outliers?)
# 2. Encontrar fuentes adicionales de variacion entre las variables (clusterizacion mas fuerte en unas variables que en otras: cepa, sexo, tratamiento...)
#    (Las variables encontradas en 2. se usaran en la design formula en los pasos siguientes)


# Unsupervised clustering analysis
vsd_t0 <- vst(dds_t0, blind = TRUE)
vsd_t1 <- vst(dds_t1, blind = TRUE)
vsd_t2 <- vst(dds_t2, blind = TRUE)
vsd_all_timebysubject <- vst(dds_all_timebysubject, blind = TRUE)
vsd_all_timetreatment <- vst(dds_all_timetreatment, blind = TRUE)

# Extraemos la matriz vst del objeto
vsd_mat_t0 <- assay(vsd_t0)
vsd_mat_t1 <- assay(vsd_t1)
vsd_mat_t2 <- assay(vsd_t2)
vsd_mat_all_timebysubject <- assay(vsd_all_timebysubject)
vsd_mat_all_timetreatment <- assay(vsd_all_timetreatment)

# Computar valores de correlacion por pares
vsd_cor_t0 <- cor(vsd_mat_t0)
vsd_cor_t1 <- cor(vsd_mat_t1)
vsd_cor_t2 <- cor(vsd_mat_t2)
vsd_cor_all_timebysubject <- cor(vsd_mat_all_timebysubject)
vsd_cor_all_timetreatment <- cor(vsd_mat_all_timetreatment)



# Plotear HEATMAP
# annotation: que factor de metadatos incluir como barras de anotacion
cor_heatmap_t0 <- pheatmap(vsd_cor_t0, 
                        annotation = dplyr::select(metadata_df_t0, treatment), 
                        main = paste0(gse_nr, " - t0"))
cor_heatmap_t1 <- pheatmap(vsd_cor_t1, 
                              annotation = dplyr::select(metadata_df_t1, treatment), 
                              main = paste0(gse_nr, " - t1"))
cor_heatmap_t2 <- pheatmap(vsd_cor_t2, 
                              annotation = dplyr::select(metadata_df_t2, treatment), 
                              main = paste0(gse_nr, " - t2"))
# En el caso del heatmap global, queremos comprobar tambien diferencia entre tejidos
cor_heatmap_all_timebysubject <- pheatmap(vsd_cor_all_timebysubject, 
                              annotation = dplyr::select(metadata_df_all, subject, time), 
                              main = paste0(gse_nr, " - Global (time and subject)"))
cor_heatmap_all_timetreatment <- pheatmap(vsd_cor_all_timetreatment, 
                                          annotation = dplyr::select(metadata_df_all, time, treatment), 
                                          main = paste0(gse_nr, " - Global (time and treatment)"))


# Guardamos heatmaps en archivos png (el global mas grande para acomodar mas muestras)
png(file.path(path_results, "correlation_heatmap_t0.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_t0
dev.off()

png(file.path(path_results, "correlation_heatmap_t1.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_t1
dev.off()

png(file.path(path_results, "correlation_heatmap_t2.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_t2
dev.off()

png(file.path(path_results, "correlation_heatmap_all_timebysubject.png"), width = 2100, height = 1750, res = 150)
cor_heatmap_all_timebysubject
dev.off()

png(file.path(path_results, "correlation_heatmap_all_timetreatment.png"), width = 2100, height = 1750, res = 150)
cor_heatmap_all_timetreatment
dev.off()

# Principal component analysis (PCA): generar plot y guardar en archivo
pca_plot_t0 <- plotPCA(vsd_t0, intgroup = "treatment") +
                  ggtitle(paste0(gse_nr, " - t0")) +
                 labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_t0.png"),
  plot = pca_plot_t0,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_t1 <- plotPCA(vsd_t1, intgroup = "treatment") +
                ggtitle(paste0(gse_nr, " - t1")) +
                labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_t1.png"),
  plot = pca_plot_t1,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_t2 <- plotPCA(vsd_t2, intgroup = "treatment") +
                 ggtitle(paste0(gse_nr, " - t2")) +
                 labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_t2.png"),
  plot = pca_plot_t2,
  width = 6,
  height = 5,
  dpi = 300
)


# A nivel global, miraremos tiempos y tratamientos
pca_plot_all_bytime <- plotPCA(vsd_all_timetreatment, intgroup = "time") +
  ggtitle(paste0(gse_nr, " - Global (by time)")) +
  labs(color = "time")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_all_bytime.png"),
  plot = pca_plot_all_bytime,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_all_timetreatment <- plotPCA(vsd_all_timetreatment, intgroup = c("time", "treatment")) +
  ggtitle(paste0(gse_nr, " - Global (by time and treatment)")) +
  labs(color = "time", shape = "treatment")
pca_plot_all_timetreatment
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_all_timetreatment.png"),
  plot = pca_plot_all_timetreatment,
  width = 6,
  height = 5,
  dpi = 300
)

# PCA 3-D
# Voy a comprobar el tercer componente principal, en un plot tridimensional
library(plotly)

# Función para crear un PCA 3D plot con plotly
pca_plot_3d <- function(vsd, metadata, grupo, titulo = "PCA 3D") {
  # Realizar PCA con prcomp
  pca <- prcomp(t(assay(vsd)))
  
  # Porcentaje de varianza explicada
  var_exp <- round(100 * summary(pca)$importance[2, 1:3], 1)
  
  # Extraer las tres primeras PCs
  pca_df <- as.data.frame(pca$x[, 1:3])
  pca_df[[grupo]] <- metadata[[grupo]]
  
  # Crear gráfico interactivo con plotly
  plot_ly(pca_df, 
          x = ~PC1, y = ~PC2, z = ~PC3, 
          color = ~get(grupo), 
          colors = "Set1",
          type = "scatter3d", 
          mode = "markers") %>%
    layout(title = titulo,
           scene = list(
             xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
             yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
             zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
           ))
}

# Ploteamos ahora para cada tejido
pca_plot_3d(vsd_t0, metadata_df_t0, "treatment", titulo = "PCA 3D - t0")
pca_plot_3d(vsd_t1, metadata_df_t1, "treatment", titulo = "PCA 3D - t1")
pca_plot_3d(vsd_t2, metadata_df_t2, "treatment", titulo = "PCA 3D - t2")
pca_plot_3d(vsd_all_timebysubject, metadata_df_all, "time", titulo = "PCA 3D - Global (by time)")


# Una nueva funcion para PCA en 3D, incluyendo dos variables (para color y forma)
pca_plot_3d <- function(vsd, metadata, color_var, shape_var = NULL, titulo = "PCA 3D") {
  # Realizar PCA con prcomp
  pca <- prcomp(t(assay(vsd)))
  
  # Porcentaje de varianza explicada
  var_exp <- round(100 * summary(pca)$importance[2, 1:3], 1)
  
  # Crear data.frame con los tres primeros PCs
  pca_df <- as.data.frame(pca$x[, 1:3])
  pca_df[[color_var]] <- metadata[[color_var]]
  pca_df[[shape_var]] <- as.factor(metadata[[shape_var]])
  
  # Crear gráfico 3D
  plot_ly(pca_df,
          x = ~PC1, y = ~PC2, z = ~PC3,
          color = ~get(color_var),
          symbol = ~get(shape_var),
          colors = "Set1",
          symbols = c("circle", "square", "diamond", "x"),  # personaliza si quieres más
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 5)) %>%
    layout(title = titulo,
           legend = list(itemsizing = "constant"),
           scene = list(
             xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
             yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
             zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
           ))
}

pca_plot_3d(vsd_all_timetreatment, metadata_df_all, color_var = "time", shape_var = "treatment",
            titulo = "PCA 3D - Global (by time and treatment)")




#### ---- ANALISIS DE EXPRESION DIFERENCIAL (DEA) ---- ####

# MODELAR RAW COUNTS PARA CADA GEN

# Contamos con varios objetos DESeq2
# dds_t0, dds_t1, dds_t2, dds_all_timebysubject, dds_all_timetreatment

# Aplicamos DESeq al DESeqDataSet
dds_t0 <- DESeq(dds_t0)
dds_t1 <- DESeq(dds_t1)
dds_t2 <- DESeq(dds_t2)
dds_all_timebysubject <- DESeq(dds_all_timebysubject)
dds_all_timetreatment <- DESeq(dds_all_timetreatment)


# PLOT DE VARIANZA VS MEDIA
# Calcular media y varianza para cada gen (cada fila)
# apply(un dataframe, 1 por filas o 2 por columnas, funcion)

# Para t0
mean_counts_t0 <- apply(exp_df_t0, 1, mean)
variance_counts_t0 <- apply(exp_df_t0, 1, var)
# Crear dataframe con media y varianza de cada gen
meanvar_df_t0 <- data.frame(mean_counts_t0, variance_counts_t0)
# Ploteamos
meanvar_plot_t0 <- ggplot(meanvar_df_t0) +
                  geom_point(aes(x=mean_counts_t0, y=variance_counts_t0)) +
                  scale_y_log10() +
                  scale_x_log10() +
                  xlab("Mean counts per gene") +
                  ylab("Variance per gene") +
                  ggtitle("Mean vs Variance – T0 Basal")
meanvar_plot_t0
ggplot2::ggsave(
  filename = file.path(path_results, "meanvar_plot_t0.png"),
  plot = meanvar_plot_t0,
  width = 6,
  height = 5,
  dpi = 300
)

# Para t1
mean_counts_t1 <- apply(exp_df_t1, 1, mean)
variance_counts_t1 <- apply(exp_df_t1, 1, var)
# Crear dataframe con media y varianza de cada gen
meanvar_df_t1 <- data.frame(mean_counts_t1, variance_counts_t1)
# Ploteamos
meanvar_plot_t1 <- ggplot(meanvar_df_t1) +
  geom_point(aes(x=mean_counts_t1, y=variance_counts_t1)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  ggtitle("Mean vs Variance – T1 Post-supplementation")
meanvar_plot_t1
ggplot2::ggsave(
  filename = file.path(path_results, "meanvar_plot_t1.png"),
  plot = meanvar_plot_t1,
  width = 6,
  height = 5,
  dpi = 300
)

# Para t2
mean_counts_t2 <- apply(exp_df_t2, 1, mean)
variance_counts_t2 <- apply(exp_df_t2, 1, var)
# Crear dataframe con media y varianza de cada gen
meanvar_df_t2 <- data.frame(mean_counts_t2, variance_counts_t2)
# Ploteamos
meanvar_plot_t2 <- ggplot(meanvar_df_t2) +
  geom_point(aes(x=mean_counts_t2, y=variance_counts_t2)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  ggtitle("Mean vs Variance – T2 Post-LPS-treatment")
meanvar_plot_t2
ggplot2::ggsave(
  filename = file.path(path_results, "meanvar_plot_t2.png"),
  plot = meanvar_plot_t2,
  width = 6,
  height = 5,
  dpi = 300
)


# Para global
mean_counts_all <- apply(exp_df_all, 1, mean)
variance_counts_all <- apply(exp_df_all, 1, var)
# Crear dataframe con media y varianza de cada gen
meanvar_df_all <- data.frame(mean_counts_all, variance_counts_all)
# Ploteamos
meanvar_plot_all <- ggplot(meanvar_df_all) +
  geom_point(aes(x=mean_counts_all, y=variance_counts_all)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  ggtitle("Mean vs Variance – Global")
meanvar_plot_all
ggplot2::ggsave(
  filename = file.path(path_results, "meanvar_plot_all.png"),
  plot = meanvar_plot_all,
  width = 6,
  height = 5,
  dpi = 300
)

# Vemos que la varianza aumenta con la media (esperable en datos RNAseq)
# La varianza tambien tiene rango mas amplio a valores bajos (esperable en datos RNAseq). La dispersion es este rango.

# Como DESeq2 modela dispersion:  Var = mean + disp * mean^2 
#    Mayor varianza -> mayor dispersion.  Mayor media -> menor dispersion.


# PLOT DE ESTIMACIONES DE DISPERSION

# Para t0
png(filename = file.path(path_results, "dispersion_plot_t0.png"), width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds_t0, main = "Dispersion estimates – T0 Basal")
dev.off()

# Para t1
png(filename = file.path(path_results, "dispersion_plot_t1.png"), width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds_t1, main = "Dispersion estimates – T1 Post-supplementation")
dev.off()

# Para t2
png(filename = file.path(path_results, "dispersion_plot_t2.png"), width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds_t2, main = "Dispersion estimates – T2 Post-LPS-treatment")
dev.off()

# Para global
png(filename = file.path(path_results, "dispersion_plot_all_timetreatment.png"), width = 6, height = 5, units = "in", res = 300)
plotDispEsts(dds_all_timetreatment, main = "Dispersion estimates – Global")
dev.off()

# Vemos que la linea de dispersion decrece con la media.
# IMPORTANTE: genes con media de expresion baja pueden dar falsos positivos!


# CALCULO DE RESULTADOS Y REDUCCION DE LOG2 FOLD CHANGES

# Calculo de resultados por defecto: test de Wald (tratamiento vs placebo, si hemos hecho relevel bien)
results_t0 <- results(dds_t0, alpha = 0.05)
results_t0
results_t1 <- results(dds_t1, alpha = 0.05)
results_t2 <- results(dds_t2, alpha = 0.05)
results_all_timebysubject <- results(dds_all_timebysubject, alpha = 0.05)
results_all_timebysubject
results_all_timetreatment <- results(dds_all_timetreatment, alpha = 0.05)
results_all_timetreatment


# Calculo de resultados con log2FC shrinkage, para conseguir log2FC mas precisos
# Version con "apeglm", menos sesgado
resultsNames(dds_t0)
results_t0_lfcshrink <- lfcShrink(dds_t0, 
                       coef = "treatment_Omega3_vs_Placebo",
                       type = "apeglm")
resultsNames(dds_t1)
results_t1_lfcshrink <- lfcShrink(dds_t1, 
                                  coef = "treatment_Omega3_vs_Placebo",
                                  type = "apeglm")
resultsNames(dds_t2)
results_t2_lfcshrink <- lfcShrink(dds_t2, 
                                  coef = "treatment_Omega3_vs_Placebo",
                                  type = "apeglm")
#resultsNames(dds_all_timebysubject)
#results_all_timebysubject_lfcshrink <- lfcShrink(dds_all_timebysubject, 
#                                  coef = "treatment_Omega3_vs_Placebo",
#                                  type = "apeglm")
#resultsNames(dds_all_timetreatment)
#results_all_timetreatment_lfcshrink <- lfcShrink(dds_all_timetreatment, 
#                                  coef = "treatment_Omega3_vs_Placebo",
#                                  type = "apeglm")

# EXPLORACION DE RESULTADOS

# Descripcion de las columnas de la tabla de resultados. as.data.frame() permite ver texto sin acortar
as.data.frame(mcols(results_t0_lfcshrink))  
as.data.frame(mcols(results_t1_lfcshrink))  
as.data.frame(mcols(results_t2_lfcshrink))  


# Observamos los resultados e identificamos los genes expresados diferencialmente (por su p-value ajustado). 
# t0
head(results_t0, n=10)
summary(results_t0)
head(results_t0_lfcshrink, n=10)
summary(results_t0_lfcshrink)
# t1
head(results_t1, n=10)
summary(results_t1)
head(results_t1_lfcshrink, n=10)
summary(results_t1_lfcshrink)
# t2
head(results_t2, n=10)
summary(results_t2)
head(results_t2_lfcshrink, n=10)
summary(results_t2_lfcshrink)


# Podemos anotar genes ahora, si no se ha hecho ya (en este caso, ya esta hecho)

# Extraemos y ordenamos los genes significativos (p-value ajustado < 0.05)
results_t2_df <- data.frame(results_t2_lfcshrink)
results_t2_df_sig <- subset(results_t2_df, padj < 0.05)
results_t2_df_sig <- results_t2_df_sig %>%
  arrange(padj)
View(results_t2_df_sig)



#### ---- VISUALIZACION DE RESULTADOS ---- ####

# EXPRESSION HEATMAP
# Ploteamos los conteos normalizados de los genes significativos

# Subset de los conteos normalizados de los genes significativos
exp_df_t2_normalizado_sig <- exp_df_t2_normalizado[rownames(results_t2_df_sig), ]
str(exp_df_t2_normalizado_sig)

# Elegir paleta de colores de la libreria RColorBrewer
# Generar gradiente de colores: azul (bajo) → blanco (0) → rojo (alto)
# breaks para truncar valores extremos de color
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-2, 2, length.out = 101)  # más resolución en torno a 0

# Expression heatmap
expr_heatmap_t2 <- pheatmap(exp_df_t2_normalizado_sig,
                   color = heat_colors,
                   breaks = breaks,
                   cluster_rows = T,
                   show_rownames = T, 
                   annotation = dplyr::select(metadata_df_t2, treatment),
                   scale = "row",
                   display_numbers = TRUE,
                   main = "Mapa de calor de expresión génica\nT2: Post-suplementación y post-tratamiento con LPS"
                  )

# Guardamos heatmaps en archivo png
png(file.path(path_results, "expression_heatmap_t2.png"), width = 1200, height = 1000, res = 150)
expr_heatmap_t2
dev.off()


# VOLCANO PLOT
# FC respecto a los p-values ajustados

# Obtenemos vector logico con los p-values ajustados < 0.05
results_t2_df <- results_t2_df %>%
  mutate(threshold = padj < 0.05) 
# %>% rownames_to_column(var = "ensgene")

# Volcano plot
volcano_plot_t2 <- ggplot(results_t2_df) +
  geom_point(aes(x=log2FoldChange, y = -log10(padj),
                 color = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
volcano_plot_t2
ggplot2::ggsave(
  filename = file.path(path_results, "volcano_plot_t2.png"),
  plot = volcano_plot_t2,
  width = 6,
  height = 5,
  dpi = 300
)


# EXPRESSION PLOT

# Seleccionamos los genes más significativos
top_genes_t2 <- data.frame(exp_df_t2_normalizado_sig)
top_genes_t2$genes <- rownames(top_genes_t2)

# Convertir las muestras en factores de una sola columna
top_genes_t2 <- gather(top_genes_t2, key = "samplename", value = "normalized_counts", 1:14)

# Unir la metadata para colorear por grupo de muestra
top_genes_t2 <- inner_join(top_genes_t2,
                     rownames_to_column(metadata_df_t2, var="samplename"),
                     by="samplename")

# Ploteamos
exp_plot_t2 <- ggplot(top_genes_t2) +
              geom_point(aes(x = genes, y = normalized_counts, color = treatment)) +
              scale_y_log10() +
              xlab("Genes") +
              ylab("Normalized Counts") +
              ggtitle("Top Significant DE Genes\nT2 Post-suplemmentation Post-LPS-treatment") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              theme(plot.title = element_text(hjust = 0.5))
exp_plot_t2
ggplot2::ggsave(
  filename = file.path(path_results, "expression_plot_t2.png"),
  plot = exp_plot_t2,
  width = 6,
  height = 5,
  dpi = 300
)

# Podemos plotear tambien por nombre de gen (mirar esto)





#### ---- RESUMEN DE CODIGO DE DATACAMP ---- ####

# Let's run through the DESeq2 workflow using the full dataset 
# with both wildtype and smoc2 overexpression samples included. 
# We have loaded the DESeq2 and dplyr libraries and read in the metadata file, 
# all_metadata and the raw counts file, all_rawcounts for you.

# Instructions
# 1. Check that the samples are in the same order in both all_rawcounts and all_metadata using the rownames(), colnames(), all(), and %in% operator.
# 2. Create the DESeq2 object using the appropriate design, testing for the effect of condition while controlling for genotype.
# 3. Create the DESeq2 object using the appropriate design, controlling for genotype and condition individually, but test for genotype:condition.

# Check that all of the samples are in the same order in the metadata and count data
all(rownames(all_metadata) %in% colnames(all_rawcounts))

# DESeq object to test for the effect of fibrosis regardless of genotype
dds_all <- DESeqDataSetFromMatrix(countData = all_rawcounts,
                                  colData = all_metadata,
                                  design = ~ genotype + condition)

# DESeq object to test for the effect of genotype on the effect of fibrosis                        
dds_complex <- DESeqDataSetFromMatrix(countData = all_rawcounts,
                                      colData = all_metadata,
                                      design = ~ genotype + condition + genotype:condition)

# Great work! We have now created our DESeq2 objects for running the differential expression analyses.


# DE analysis

# We are going to continue using the full dataset comparing the genes 
# that exhibit significant expression differences between normal and fibrosis samples 
# regardless of genotype (design: ~ genotype + condition). Therefore, we will use 
# our dds_all DESeq2 object created in the previous exercise. 
# Assume this object is created and all libraries are loaded. 
# In this exercise let's perform the unsupervised clustering analyses 
# to explore the clustering of our samples and sources of variation.

# Instructions
# Log transform the normalized counts inside the dds_all object using the vst() function, being blind to sample group information.
# Create the correlation heatmap of the correlation values of the log normalized counts using the pheatmap() function. Include annotation bars for genotype and condition.
# Plot the PCA with the plotPCA() function using vsd_all. Color the plot by condition.
# Plot the PCA with the plotPCA() function using vsd_all. Color the plot by genotype.

# Log transform counts for QC
vsd_all <- vst(dds_all, blind = TRUE)

# Create heatmap of sample correlation values
vsd_all %>% 
  assay() %>%
  cor() %>%
  pheatmap(annotation = select(all_metadata, c("genotype", "condition")))

# Create the PCA plot for PC1 and PC2 and color by condition       
plotPCA(vsd_all, intgroup  = "condition")

# Create the PCA plot for PC1 and PC2 and color by genotype       
plotPCA(vsd_all, intgroup  = "genotype")

# Awesome work! Now that we have explored the data and identified 
# sources of variation and clustering of our samples we are ready for DE analysis.


# DE analysis results

# After exploring the PCA and correlation heatmap, we found good clustering 
# of our samples on PC1, which seemed to represent the variation in the data 
# due to fibrosis, and PC2, which appeared to represent variation in the data 
# due to smoc2 overexpression. We did not find additional sources of variation 
# in the data, nor any outliers to remove. Therefore, we can proceed by 
# running DESeq2, DE testing, and shrinking the fold changes. We performed 
# these steps for you to generate the final results, res_all.

# In this exercise, we'll want to subset the significant genes from the results 
# and output the top 10 DE genes by adjusted p-value.

# Instructions
# Use the subset() function to extract those values with an adjusted p-value less than 0.05. Save the subset as a data frame named smoc2_sig by using the data.frame() function and turning the row names to a column named geneID using the rownames_to_column() function.
# Order the significant results by adjusted p-values using the arrange() function, select the columns with Ensembl gene ID and adjusted p-values, and output the top significant genes using head().

# Select significant genese with padj < 0.05
smoc2_sig <- subset(res_all, padj < 0.05) %>%
  data.frame() %>%
  rownames_to_column(var = "geneID")

# Extract the top 6 genes with padj values
smoc2_sig %>%
  arrange(padj) %>%
  select(geneID, padj) %>%
  head()

# Great work! You have identified the differentially expressed genes. Now it would 
# be time to validate any of the interesting results! Also, remember that we can use 
# the annotables package to convert the Ensembl IDs to the more recognizable gene symbols.
