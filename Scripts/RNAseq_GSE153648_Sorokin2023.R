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
gse_nr <- "GSE153648"
gse_description <- "GSE153648_Sorokin2023"

# Marcamos la ruta y guardamos nombres de archivo en una lista
path_ubuntu <- paste0("//wsl.localhost/Ubuntu/home/victor/UNIR/TFM/Datasets/", gse_description, "/")
setwd(path_ubuntu)
file_names <- list.files(path = path_ubuntu, pattern = "\\.xlsx$")
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

# Cargamos readxl para leer archivo excel
library(readxl)
exp_df_all <- read_excel(file_path, sheet=1)

# Extraemos dos primeras columnas, con nombres y codigos de genes, en un dataframe aparte
genes <- exp_df_all[,c("Name", "Identifier")]

# Eliminamos las columnas de genes y 
exp_df_all <- exp_df_all %>% dplyr::select(-c("Name", "Identifier"))
rownames(exp_df_all) <- genes$Name

# Exploramos rapidamente las primeras filas y la estructura del dataframe de expresion
head(exp_df_all)
str(exp_df_all)



# CREAMOS DATAFRAME DE METADATOS

# Acceso a metadatos usando getGEO
gse <- getGEO(gse_nr, GSEMatrix = TRUE)
names(gse)
names(pData(gse[[1]]))
metadata_df_all <- pData(phenoData(gse[[1]]))

# Simplificamos nombre de la columna tissue
metadata_df_all$tissue <- metadata_df_all[["tissue:ch1"]]
metadata_df_all[["tissue:ch1"]] <- NULL
# Simplificamos la columna tratamiento (nombre, y los espacios en blanco por -)
metadata_df_all$treatment <- metadata_df_all[["treatment:ch1"]]
metadata_df_all[["treatment:ch1"]] <- NULL
metadata_df_all$treatment <- gsub(" ", "-", metadata_df_all$treatment)



# AJUSTE DE DATOS Y METADATOS
# (Probablemente algunas lineas aqui son redundantes)

# Extraemos los nombres C## a una nueva columna
metadata_df_all <- metadata_df_all %>%
  mutate(C_muestra = str_extract(title, "(?<=, ).*"))

# Columna de nombres de muestra
sample_name <- colnames(exp_df_all)
sample_name

# Simplificamos los nombres de muestras y sustituimos nombres de muestra en dataframe de conteos
sample_name <- str_extract(sample_name, "(?<=-).*(?= -)")
colnames(exp_df_all) <- sample_name

# Guardamos los nombres GSM de metadatos en el mismo orden que los nombres C de muestras
sample_name_new <- metadata_df_all$geo_accession[match(sample_name, metadata_df_all$C_muestra)]

# Con funcion match() reordenamos muestras en matriz de conteos para que coincida con el orden en metadatos
reorder_idx <- match(metadata_df_all$C_muestra, colnames(exp_df_all))
exp_df_all <- exp_df_all[ , reorder_idx]

# Cambiamos ahora los nombres C por los nombres GSM en el orden correcto
colnames(exp_df_all) <- rownames(metadata_df_all)
rownames(exp_df_all) <- genes$Name

# Nos aseguramos de que los nombres de muestras estan en el mismo orden como columnas de datos y filas de metadatos
print(paste0("Coinciden los nombres de todas las muestras: ", all(row.names(metadata_df_all) == colnames(exp_df_all))))
row.names(metadata_df_all) == colnames(exp_df_all)


# DIVISION DEL DATASET EN SUB-DATASETS

# Numeros de muestra para cada tejido
samples_aorta <- paste0("GSM", seq(4648593, 4648606))
samples_skin <- paste0("GSM", seq(4648607, 4648620))
samples_liver <- paste0("GSM", seq(4648621, 4648634))

# Metadata de cada tejido
metadata_df_aorta <- metadata_df_all[samples_aorta, ]
metadata_df_skin <- metadata_df_all[samples_skin, ]
metadata_df_liver <- metadata_df_all[samples_liver, ]

# Matriz de conteo de cada tejido
exp_df_aorta <- exp_df_all[,samples_aorta]
exp_df_skin <- exp_df_all[,samples_skin]
exp_df_liver <- exp_df_all[,samples_liver]


# FACTORES

# Creamos factores para todas las variables de estudio
treatment_all <- factor(metadata_df_all$treatment) 
names(treatment_all) <- names(exp_df_all) # Mantenemos nombres de muestras en el factor grupos
treatment_all <- relevel(treatment_all, ref = "Omega-3-Deficient") # Nos aseguramos de que Placebo (Control) es el factor de referencia
# Esto ayuda a no cometer errores a la hora de hacer comparaciones de grupos

# Dividimos factores de los tres tejidos
treatment_aorta <- treatment_all[samples_aorta]
treatment_skin <- treatment_all[samples_skin]
treatment_liver <- treatment_all[samples_liver]


### --- OBJETO DESEQDATASET

# OBJETOS DESEQDATASET
dds_aorta <- DESeqDataSetFromMatrix(countData =  exp_df_aorta,
                                  colData =  metadata_df_aorta,
                                  design = ~ treatment)
dds_skin <- DESeqDataSetFromMatrix(countData =  exp_df_skin,
                                    colData =  metadata_df_skin,
                                    design = ~ treatment)
dds_liver <- DESeqDataSetFromMatrix(countData =  exp_df_liver,
                                    colData =  metadata_df_liver,
                                    design = ~ treatment)
dds_all <- DESeqDataSetFromMatrix(countData =  exp_df_all,
                                    colData =  metadata_df_all,
                                    design = ~ treatment + tissue)


# NORMALIZACION DE LOS CONTEOS

# Determinar los size factors a usar en la normalizacion
dds_aorta <- estimateSizeFactors(dds_aorta)
dds_skin <- estimateSizeFactors(dds_skin)
dds_liver <- estimateSizeFactors(dds_liver)
dds_all <- estimateSizeFactors(dds_all)

# Extraer el dataframe de conteos normalizados
exp_df_aorta_normalizado <- counts(dds_aorta, normalized=TRUE)
exp_df_skin_normalizado <- counts(dds_skin, normalized=TRUE)
exp_df_liver_normalizado <- counts(dds_liver, normalized=TRUE)
exp_df_all_normalizado <- counts(dds_all, normalized=TRUE)



#### ---- EXPLORACION DE LOS DATOS ---- ####

# CLUSTERIZACION NO SUPERVISADA (heatmap y PCA)
# 1. Hacer un analisis de calidad de las muestras (se agrupan muestras similares entre si como esperamos? hay outliers?)
# 2. Encontrar fuentes adicionales de variacion entre las variables (clusterizacion mas fuerte en unas variables que en otras: cepa, sexo, tratamiento...)
#    (Las variables encontradas en 2. se usaran en la design formula en los pasos siguientes)


# Unsupervised clustering analysis
vsd_aorta <- vst(dds_aorta, blind = TRUE)
vsd_skin <- vst(dds_skin, blind = TRUE)
vsd_liver <- vst(dds_liver, blind = TRUE)
vsd_all <- vst(dds_all, blind = TRUE)

# Extraemos la matriz vst del objeto
vsd_mat_aorta <- assay(vsd_aorta)
vsd_mat_skin <- assay(vsd_skin)
vsd_mat_liver <- assay(vsd_liver)
vsd_mat_all <- assay(vsd_all)

# Computar valores de correlacion por pares
vsd_cor_aorta <- cor(vsd_mat_aorta)
vsd_cor_skin <- cor(vsd_mat_skin)
vsd_cor_liver <- cor(vsd_mat_liver)
vsd_cor_all <- cor(vsd_mat_all)



# Plotear HEATMAP
# annotation: que factor de metadatos incluir como barras de anotacion
cor_heatmap_aorta <- pheatmap(vsd_cor_aorta, 
                        annotation = dplyr::select(metadata_df_aorta, treatment), 
                        main = paste0(gse_nr, " - Aorta"))
cor_heatmap_skin <- pheatmap(vsd_cor_skin, 
                              annotation = dplyr::select(metadata_df_skin, treatment), 
                              main = paste0(gse_nr, " - Skin"))
cor_heatmap_liver <- pheatmap(vsd_cor_liver, 
                              annotation = dplyr::select(metadata_df_liver, treatment), 
                              main = paste0(gse_nr, " - Liver"))
# En el caso del heatmap global, queremos comprobar tambien diferencia entre tejidos
cor_heatmap_all <- pheatmap(vsd_cor_all, 
                              annotation = dplyr::select(metadata_df_all, tissue, treatment), 
                              main = paste0(gse_nr, " - Global"))


# Guardamos heatmaps en archivos png (el global mas grande para acomodar mas muestras)
png(file.path(path_results, "correlation_heatmap_aorta.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_aorta
dev.off()

png(file.path(path_results, "correlation_heatmap_skin.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_skin
dev.off()

png(file.path(path_results, "correlation_heatmap_liver.png"), width = 1200, height = 1000, res = 150)
cor_heatmap_liver
dev.off()

png(file.path(path_results, "correlation_heatmap_all.png"), width = 2100, height = 1750, res = 150)
cor_heatmap_all
dev.off()


# Principal component analysis (PCA): generar plot y guardar en archivo
pca_plot_aorta <- plotPCA(vsd_aorta, intgroup = "treatment") +
                  ggtitle(paste0(gse_nr, " - Aorta")) +
                 labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_aorta.png"),
  plot = pca_plot_aorta,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_skin <- plotPCA(vsd_skin, intgroup = "treatment") +
                ggtitle(paste0(gse_nr, " - Skin")) +
                labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_skin.png"),
  plot = pca_plot_skin,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_liver <- plotPCA(vsd_liver, intgroup = "treatment") +
                 ggtitle(paste0(gse_nr, " - Liver")) +
                 labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_liver.png"),
  plot = pca_plot_liver,
  width = 6,
  height = 5,
  dpi = 300
)


# A nivel global, miraremos tejidos y tratamientos
pca_plot_all_bytissue <- plotPCA(vsd_all, intgroup = "tissue") +
  ggtitle(paste0(gse_nr, " - Global (by tissue)")) +
  labs(color = "tissue")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_all_bytissue.png"),
  plot = pca_plot_all_bytissue,
  width = 6,
  height = 5,
  dpi = 300
)

pca_plot_all_bytreatment <- plotPCA(vsd_all, intgroup = "treatment") +
  ggtitle(paste0(gse_nr, " - Global (by treatment)")) +
  labs(color = "treatment")
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot_all_bytreatment.png"),
  plot = pca_plot_all_bytreatment,
  width = 6,
  height = 5,
  dpi = 300
)

# Voy a comprobar el tercer componente principal, en un plot tridimensional
library(plotly)

# Función para crear un PCA 3D plot con plotly
pca_plot_3d <- function(vsd, metadata, grupo, titulo = "PCA 3D") {
  # Realizar PCA con prcomp
  pca <- prcomp(t(assay(vsd)))
  
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
           scene = list(xaxis = list(title = "PC1"),
                        yaxis = list(title = "PC2"),
                        zaxis = list(title = "PC3")))
}

# Ploteamos ahora para cada tejido
pca_plot_3d()

# Control de calidad no superado.
# FIN DEL CODIGO