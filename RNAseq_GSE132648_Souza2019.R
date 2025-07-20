rm(list=ls())


# CARGAMOS LIBRERIAS

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)

# Para informacion detallada de DESeq2, ejecutar esta linea:
#vignette("DESeq2")


#### ---- CARGA DE DATOS ---- ####

# NOTAS: 
# DATOS: Es importante usar CONTEOS NO NORMALIZADOS como input para DESeq2 (no normalizados, no transformados)
# DATOS vs METADATOS: los nombres de muestras deben encontrarse como nombres de columnas en los DATOS y como nombres de filas en los METADATOS, y en el mismo orden.

# CARGAMOS RUTA Y NOMBRES DE ARCHIVO

# Indicar codigo de experimento y descripcion (p.ej. con autor). Para rutas de archivo
gse_nr <- "GSE132648"
gse_description <- "GSE132648_Souza2019"

# Marcamos la ruta y guardamos nombres de archivo en una lista
path_ubuntu <- paste0("//wsl.localhost/Ubuntu/home/victor/UNIR/TFM/Datasets/", gse_description, "/")
setwd(path_ubuntu)
file_names <- list.files(path = path_ubuntu, pattern = "\\.txt$")
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
exp_df_all <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

# Convertimos la primera columna en rownames pero asegurando que sean únicos
exp_df_all <- exp_df_all[!duplicated(exp_df_all$X), ]
row.names(exp_df_all) <- exp_df_all$X

# Eliminamos la primera columna (ya está como rownames)
exp_df_all <- exp_df_all[, -1]


# Exploramos rapidamente las primeras filas y la estructura del dataframe de expresion
head(exp_df_all)
str(exp_df_all)



# CREAMOS DATAFRAME DE METADATOS

# Columna de nombres de muestra
sample_name <- colnames(exp_df_all)
sample_name


# Abrimos archivo de metadatos
metadata_path <- paste0(path_ubuntu, "GSE132648_metadata.txt")
metadata_df_all <- read.delim(metadata_path, header = TRUE, stringsAsFactors = FALSE)

# Ahora puedes ver o trabajar con los metadatos
head(metadata_df_all)

# Mantenemos solo muestras con quality ok
metadata_df_all <- metadata_df_all[metadata_df_all$quality == "ok", ]

# Cambiamos muestras de los metadatos a nombres de filas
row.names(metadata_df_all) <- metadata_df_all$sample
metadata_df_all <- metadata_df_all[,-1]


# Nos aseguramos de que los nombres de muestras estan en el mismo orden como columnas de datos y filas de metadatos
print(paste0("Coinciden los nombres de todas las muestras: ", all(row.names(metadata_df_all) == colnames(exp_df_all))))
row.names(metadata_df_all) == colnames(exp_df_all)

# Si no coinciden los nombres podemos usar este código:
# Usamos la funcion match() para reordenar las columnas de los conteos
indices_reorden <- match(rownames(metadata_df_all), colnames(exp_df_all))
# Reordenamos las columnas de la matriz de conteos
exp_df_all <- exp_df_all[ , indices_reorden]


# Creamos factores para todas las variables de estudio
groups <- factor(metadata_df_all$group) 
names(groups) <- names(exp_df_all) # Mantenemos nombres de muestras en el factor grupos
groups <- relevel(groups, ref = "Placebo") # Nos aseguramos de que Placebo (Control) es el factor de referencia
# Esto ayuda a no cometer errores a la hora de hacer comparaciones de grupos

batches <- factor(metadata_df_all$batch) # Factor batch
subjects <- factor(metadata_df_all$subject) # Sujetos como factor para el analisis en pares 



# OBJETO DESEQDATASET
dds_all <- DESeqDataSetFromMatrix(countData =  exp_df_all,
                                  colData =  metadata_df_all,
                                  design = ~ group)

# NORMALIZACION DE LOS CONTEOS

# Determinar los size factors a usar en la normalizacion
dds_all <- estimateSizeFactors(dds_all)

# Extraer el dataframe de conteos normalizados
exp_df_normalizado <- counts(dds_all, normalized=TRUE)



#### ---- EXPLORACION DE LOS DATOS ---- ####

# CLUSTERIZACION NO SUPERVISADA (heatmap y PCA)
# 1. Hacer un analisis de calidad de las muestras (se agrupan muestras similares entre si como esperamos? hay outliers?)
# 2. Encontrar fuentes adicionales de variacion entre las variables (clusterizacion mas fuerte en unas variables que en otras: cepa, sexo, tratamiento...)
#    (Las variables encontradas en 2. se usaran en la design formula en los pasos siguientes)


# Unsupervised clustering analysis
vsd_all <- vst(dds_all, blind = TRUE)

# Extraemos la matriz vst del objeto
vsd_mat_all <- assay(vsd_all)

# Computar valores de correlacion por pares
vsd_cor_all <- cor(vsd_mat_all)
View(vsd_cor_all)

# Plotear HEATMAP
# annotation: que factor de metadatos incluir como barras de anotacion
cor_heatmap <- pheatmap(vsd_cor_all, 
                        annotation = dplyr::select(metadata_df_all, group), 
                        main = gse_nr)

# Guardamos heatmap en un archivo png
png(file.path(path_results, "correlation_heatmap.png"), width = 1200, height = 1000, res = 150)
cor_heatmap
dev.off()


# Principal component analysis (PCA): generar plot y guardar en archivo
pca_plot <- plotPCA(vsd_all, intgroup = "group") +
            ggtitle(gse_nr)
ggplot2::ggsave(
  filename = file.path(path_results, "pca_plot.png"),
  plot = pca_plot,
  width = 6,
  height = 5,
  dpi = 300
)


# Control de calidad no superado.
# FIN DEL CODIGO