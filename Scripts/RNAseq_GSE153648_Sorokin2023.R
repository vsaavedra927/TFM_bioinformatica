# Resetear
rm(list = ls())           # limpiar workspace


# Establecer ruta de la carpeta
path_carpetaGithub <- "C:/Users/victo/OneDrive/Escritorio/Github"

# Establecer valores de interes
umbral_fdr <- 0.05
umbral_lfc <- 0.585





# Cargar funciones y paquetes del archivo de funciones
path_funciones <- paste0(path_carpetaGithub, "/Scripts")
source(paste0(path_funciones, "/RNAseq_funciones.R"))
load_pkgs()




#/////////////////////////////////////
#---- GESTION SUBSETS Y METADATOS ----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



#### ---- CARGA DE DATOS ---- ####

# NOTAS: 
# DATOS: Es importante usar CONTEOS NO NORMALIZADOS como input para DESeq2 (no normalizados, no transformados)
# DATOS vs METADATOS: los nombres de muestras deben encontrarse como nombres de columnas en los DATOS y como nombres de filas en los METADATOS, y en el mismo orden.

# CARGAMOS RUTA Y NOMBRES DE ARCHIVO

# Indicar codigo de experimento y descripcion (p.ej. con autor). Para rutas de archivo
gse_nr <- "GSE153648"
gse_description <- "GSE153648_Sorokin2023"
especie <- "Mus musculus"

# Marcamos la ruta y guardamos nombres de archivo en una lista
path_dataset <- paste0(path_carpetaGithub, "/Datasets/", gse_description, "/")
file_names <- list.files(path = path_dataset, pattern = "\\.xlsx$")
file_names
data_file <- file_names[1]

# Marcamos la ruta para guardar resultados. Crearla si no existe
path_results <- paste0(path_carpetaGithub, "/Results/Results_", gse_description, "/")
if (!dir.exists(path_results)) {
  dir.create(path_results, recursive = TRUE)
}



# CREAMOS DATAFRAME DE EXPRESION

# Existe un solo archivo para abrir en un dataframe

# Construimos ruta del archivo
file_path <- paste0(path_dataset, data_file)

# Cargamos readxl para leer archivo excel
library(readxl)
# Cargamos archivo, y que todos los valores numéricos se lean correctamente como numeric
exp_df_all <- read_excel(file_path, sheet = 1, col_types = "text") %>%
  type_convert()


# Extraemos dos primeras columnas, con nombres y codigos de genes, en un dataframe aparte
id_genes <- exp_df_all[,c("Name", "Identifier")]

# Eliminamos las columnas de genes y 
exp_df_all <- exp_df_all %>% dplyr::select(-c("Name", "Identifier"))
rownames(exp_df_all) <- id_genes$Name

# Exploramos rapidamente las primeras filas y la estructura del dataframe de expresion
head(exp_df_all)
str(exp_df_all)




# CREAMOS DATAFRAME DE METADATOS

# Acceso a metadatos usando getGEO
gse <- getGEO(gse_nr, GSEMatrix = TRUE)
names(gse)
names(pData(gse[[1]]))
metadata_df_all <- pData(phenoData(gse[[1]]))

# Simplificamos nombre de la columna tejido
metadata_df_all$tejido <- metadata_df_all[["tissue:ch1"]]
metadata_df_all[["tissue:ch1"]] <- NULL
# Simplificamos la columna tratamiento (nombre y valores)
metadata_df_all$tratamiento <- metadata_df_all[["treatment:ch1"]]
metadata_df_all[["treatment:ch1"]] <- NULL
metadata_df_all <- metadata_df_all %>%
  mutate(
    tratamiento = case_when(
      tratamiento == "Omega-3 Deficient" ~ "Control",
      tratamiento == "EPA-suppl" ~ "EPAsup",
      tratamiento == "DHA-suppl" ~ "DHAsup",
      TRUE ~ tratamiento   # Deja otros valores como estan
    ),
    tejido = case_when(
      tejido == "aorta" ~ "Aorta",
      tejido == "skin" ~ "Piel",
      tejido == "liver" ~ "Hígado",
      TRUE ~ tejido
    )
  )



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
rownames(exp_df_all) <- id_genes$Name

# Nos aseguramos de que los nombres de muestras estan en el mismo orden como columnas de datos y filas de metadatos
print(paste0("Coinciden los nombres de todas las muestras: ", all(row.names(metadata_df_all) == colnames(exp_df_all))))
row.names(metadata_df_all) == colnames(exp_df_all)


# CONVERTIR COLUMNAS DE INTERES A FACTOR

# Convertimos variables de interes a factor y marcamos nivel de referencia
metadata_df_all$tratamiento <- factor(metadata_df_all$tratamiento) 
metadata_df_all$tratamiento <- relevel(metadata_df_all$tratamiento, ref = "Control")
levels(metadata_df_all$tratamiento)


# CONVERTIMOS DATAFRAME A MATRIZ
count_matrix_all <- as.matrix(exp_df_all)


# GENERAR SUBSETS DE DATOS EN FORMATO LISTA

# Generamos lista de propiedades del set completo
set_all <- list(nombre_subset = "all",
                gse_nr = gse_nr,
                especie = especie,
                tejido = "Global",
                tratamiento = c("Control", "DHAsup", "EPAsup"),
                muestras = rownames(metadata_df_all),
                metadata_df = metadata_df_all,
                count_matrix = count_matrix_all
)

# Crear los subsets de muestras de interes, con sus propiedades basicas
# nombre_subset: para referenciar, usado como entrada en las funciones
# gse_nr y especie: los incorporamos más abajo
# tejido: valor de la columna tejido
# tratamiento: valores de la columna tratamiento

if(TRUE){
    subset_aorta_all <- list(nombre_subset = "aorta_all",
                             tejido = "Aorta",
                             tratamiento = c("Control", "DHAsup", "EPAsup")
    )
    subset_aorta_dhavscontrol <- list(nombre_subset = "aorta_dhavscontrol",
                                      tejido = "Aorta",
                                      tratamiento = c("Control", "DHAsup")
    )
    subset_aorta_epavscontrol <- list(nombre_subset = "aorta_epavscontrol",
                                      tejido = "Aorta",
                                      tratamiento = c("Control", "EPAsup")
    )
    subset_aorta_epavsdha <- list(nombre_subset = "aorta_epavsdha",
                                  tejido = "Aorta",
                                  tratamiento = c("DHAsup", "EPAsup")
    )
    
    subset_piel_all <- list(nombre_subset = "piel_all",
                            tejido = "Piel",
                            tratamiento = c("Control", "DHAsup", "EPAsup")
    )
    subset_piel_dhavscontrol <- list(nombre_subset = "piel_dhavscontrol",
                                     tejido = "Piel",
                                     tratamiento = c("Control", "DHAsup")
    )
    subset_piel_epavscontrol <- list(nombre_subset = "piel_epavscontrol",
                                     tejido = "Piel",
                                     tratamiento = c("Control", "EPAsup")
    )
    subset_piel_epavsdha <- list(nombre_subset = "piel_epavsdha",
                                 tejido = "Piel",
                                 tratamiento = c("DHAsup", "EPAsup")
    )
    
    subset_higado_all <- list(nombre_subset = "higado_all",
                              tejido = "Hígado",
                              tratamiento = c("Control", "DHAsup", "EPAsup")
    )
    subset_higado_dhavscontrol <- list(nombre_subset = "higado_dhavscontrol",
                                       tejido = "Hígado",
                                       tratamiento = c("Control", "DHAsup")
    )
    subset_higado_epavscontrol <- list(nombre_subset = "higado_epavscontrol",
                                       tejido = "Hígado",
                                       tratamiento = c("Control", "EPAsup")
    )
    subset_higado_epavsdha <- list(nombre_subset = "higado_epavsdha",
                                   tejido = "Hígado",
                                   tratamiento = c("DHAsup", "EPAsup")
    )
}


nombres_subsets <- c("subset_aorta_all", "subset_aorta_dhavscontrol", "subset_aorta_epavscontrol", "subset_aorta_epavsdha",
                     "subset_piel_all", "subset_piel_dhavscontrol", "subset_piel_epavscontrol", "subset_piel_epavsdha", 
                     "subset_higado_all", "subset_higado_dhavscontrol", "subset_higado_epavscontrol", "subset_higado_epavsdha")

# (Debugging) Codigo para eliminar completamente los subsets 
#for (nombre_subset  in nombres_subsets){rm(list = nombre)}

# ASIGNAR VALORES DE CADA SUBSET
for (nombre_subset in nombres_subsets) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Asignar valores
  subset_obj <- asignar_gse_y_especie(subset_obj, 
                                      gse_nr = gse_nr, 
                                      especie = especie)
  subset_obj <- asignar_muestras(subset_obj, metadatos_global = metadata_df_all)
  subset_obj <- asignar_metadatos(subset_obj, metadatos_global = metadata_df_all)
  subset_obj <- asignar_matriz_conteos(subset_obj, count_matrix_global = count_matrix_all)
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}



# COLORES
colores_tejido <- setNames(c("#A65628", "#FFFF33", "#984EA3"), unique(metadata_df_all$tejido))
colores_tratamiento <- setNames(c("#999999", "#B35C00", "#FFC300"), unique(metadata_df_all$tratamiento))

# Creamos lista con las paletas de colores
set_colores = list(
  tejido = colores_tejido,
  tratamiento = colores_tratamiento
)


#-----------------------------------------------------------------
#-----------------------------------------------------------------
  
  
#/////////////////////////////////////////
#------------ ANALISIS DESEQ2 ------------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#                    Y
#//////////////////////////////////////////
#---- EXPLORACION PREVIA DE LOS DATOS -----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



# OBJETO DESEQDATASET, NORMALIZACION DE CONTEOS, EXPLORACION PREVIA DE DATOS
# Aplicar a cada subset
for (nombre_subset in nombres_subsets) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  subset_obj <- crear_DESeqDataSet(subset_obj, variable_diseno = "tratamiento")
  subset_obj <- normalizar_conteos(subset_obj)
  subset_obj <- estudiar_correlacion_datos(subset_obj, 
                                           heatmap = TRUE, pca = TRUE, 
                                           colores = set_colores)
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}



# APLICAR A SET GLOBAL (con IF para ejecutar en bloque)
# (Heatmap y PCA un poco diferentes, uso funcion adaptada)
if (TRUE){
  set_all <- crear_DESeqDataSet(set_all, variable_diseno = "tejido + tratamiento")
  set_all <- normalizar_conteos(set_all)
  # Solo correlacion
  set_all <- estudiar_correlacion_datos(set_all, 
                                           heatmap = FALSE, pca = FALSE, 
                                           colores = set_colores)
  
  # Plotear HEATMAP global
  width <- 2100
  height <- 1750
  res <- 150
  # Indicar título y ruta de archivo
  etiqueta <- paste0("Mapa de calor de correlación de muestras\n", gse_nr, " (", especie, ") - ", set_all$tejido, "\n", paste(set_all$tratamiento, collapse = ", "))
  file_name <- paste0("correlation_heatmap_all.png")
  pheatmap(set_all$vsd_cor, 
           annotation = dplyr::select(set_all$metadata_df, tejido, tratamiento), 
           annotation_colors = set_colores,
           main = etiqueta,
           filename = file.path(path_results, file_name),
           width = width/res, height = height/res   # pheatmap trabaja en pulgadas; pulgadas = pixeles / resolucion
  )
  message(paste0("Guardado heatmap de correlacion: ", "correlation_heatmap_all.png"))
  
  # PCA global de variable tratamiento
  # Generar plot
  titulo <- "Análisis de componentes principales (PCA)"
  subtitulo <- paste0(gse_nr, " (", especie, ") - ", set_all$tejido, "\n", paste(set_all$tratamiento, collapse = ", "))
  set_all$pca_plot_portratamiento <- plotPCA(set_all$vsd, intgroup = "tratamiento") +
    ggtitle(titulo, subtitle = subtitulo) +
    labs(color = "tratamiento") +
    scale_color_manual(values = set_colores$tratamiento)
  # Guardar en archivo
  file_name <- paste0("pca_plot_", set_all$nombre_subset , "_portratamiento.png")
  ggplot2::ggsave(
    filename = file.path(path_results, file_name),
    plot = set_all$pca_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
  message(paste0("Guardado PCA: ", file_name))
  # PCA global de variable tejido
  # Generar plot
  titulo <- "Análisis de componentes principales (PCA)"
  subtitulo <- paste0(gse_nr, " (", especie, ") - ", set_all$tejido, "\n", paste(set_all$tratamiento, collapse = ", "))
  set_all$pca_plot_portejido <- plotPCA(set_all$vsd, intgroup = "tejido") +
    ggtitle(titulo, subtitle = subtitulo) +
    labs(color = "tejido") +
    scale_color_manual(values = set_colores$tejido)
  # Guardar en archivo
  file_name <- paste0("pca_plot_", set_all$nombre_subset , "_portejido.png")
  ggplot2::ggsave(
    filename = file.path(path_results, file_name),
    plot = set_all$pca_plot_portejido,
    width = 6,
    height = 5,
    dpi = 300
  )
  message(paste0("Guardado PCA: ", file_name))
}


# Ploteamos ahora PCA 3D para diferenciar tejidos
etiqueta <- paste0("Análisis de componentes principales (PCA) 3D\n", gse_nr, " (", especie, ") - ", set_all$tejido, "\n", paste(set_all$tratamiento, collapse = ", "))
pca_plot_3d(set_all$vsd, set_all$metadata_df, "tratamiento", colores = set_colores$tratamiento, titulo = "PCA 3D - Global")
pca_plot_3d(set_all$vsd, set_all$metadata_df, "tejido", colores = set_colores$tejido, titulo = "PCA 3D - Global")






#//////////////////////////////////////////////////////////////
#---- ANALISIS DE EXPRESION DIFERENCIAL (DEA)  CON DESEQ2 -----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

reducir_lfc <- FALSE

# MODELAR EXPRESION DIFERENCIAL, GRAFICAR MEDIA VS VARIANZA
# Aplicar a cada subset
for (nombre_subset in nombres_subsets) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  subset_obj <- aplicar_deseq(subset_obj)
  plot_varvsmean(subset_obj)
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}





# CALCULAR RESULTADOS POR DEFECTO (SIN REDUCCION DE LOG2 FOLD CHANGES)
# NOTA: nos centraremos en los subsets que tienen un unico par de grupos a comparar

# Mostrar posible coeficientes de cada subset
for (nombre_subset in c(nombres_subsets)) {
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  cat("\nResultados para", nombre_subset, ":\n")
  print(resultsNames(subset_obj$dds))
}

# Enumerar todos los subsets que vamos a analizar sin reduccion de LFC
subsets_para_resultados <- c("subset_aorta_dhavscontrol", "subset_aorta_epavscontrol", "subset_aorta_epavsdha",
                                        "subset_piel_dhavscontrol", "subset_piel_epavscontrol", "subset_piel_epavsdha",
                                        "subset_higado_dhavscontrol", "subset_higado_epavscontrol", "subset_higado_epavsdha")

# Extraer resultados de cada subset, con el coeficiente de interes
if(TRUE){
  # aorta
  subset_aorta_dhavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_aorta_dhavscontrol, 
                                  coeficiente = "tratamiento_DHAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_aorta_epavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_aorta_epavscontrol, 
                                  coeficiente = "tratamiento_EPAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_aorta_epavsdha <- extraer_resultados_log2fcshrink(
                                  subset_aorta_epavsdha, 
                                  coeficiente = "tratamiento_EPAsup_vs_DHAsup", 
                                  control_group = "DHAsup",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  
  # piel
  subset_piel_dhavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_piel_dhavscontrol, 
                                  coeficiente = "tratamiento_DHAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_piel_epavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_piel_epavscontrol, 
                                  coeficiente = "tratamiento_EPAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_piel_epavsdha <- extraer_resultados_log2fcshrink(
                                  subset_piel_epavsdha, 
                                  coeficiente = "tratamiento_EPAsup_vs_DHAsup", 
                                  control_group = "DHAsup",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  
  # higado
  subset_higado_dhavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_higado_dhavscontrol, 
                                  coeficiente = "tratamiento_DHAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_higado_epavscontrol <- extraer_resultados_log2fcshrink(
                                  subset_higado_epavscontrol, 
                                  coeficiente = "tratamiento_EPAsup_vs_Control", 
                                  control_group = "Control",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)
  
  subset_higado_epavsdha <- extraer_resultados_log2fcshrink(
                                  subset_higado_epavsdha, 
                                  coeficiente = "tratamiento_EPAsup_vs_DHAsup", 
                                  control_group = "DHAsup",
                                  metadatos_global = metadata_df_all,
                                  reducir_lfc = reducir_lfc, 
                                  umbral_fdr = umbral_fdr,
                                  umbral_lfc = umbral_lfc)

}




#//////////////////////////////////////////////////////////////////////////////////////////////
#---------- VISUALIZACION DE RESULTADOS DEL ANALISIS DE EXPRESION DIFERENCIAL (DEA) -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# Indicar las columnas de metadatos con interes para el expression plot (si se duda, incluir mas columnas)
col_metadata_para_exp_plot <- c("title","geo_accession","organism_ch1","strain:ch1",
                                "tratamiento","tejido","C_muestra")

# VOLCANO PLOT DE SUBSETS SIN REDUCCION LFC
for (nombre_subset in c(subsets_para_resultados)) {
  print(nombre_subset)  
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  subset_obj <- visualizar_resultados(
                      subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
                      umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
                      colores = set_colores, jit_w = 0,
                      heatmap = FALSE, heatmap_labels = FALSE, volcano_plot = TRUE,
                      exp_plot = FALSE, dibujar_violin = FALSE, exp_plot_emparejado = FALSE
                    )
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}

# HEATMAPS CON ETIQUETAS DE GENES PEQUEÑAS
for (nombre_subset in c(subsets_para_resultados)) {
  print(nombre_subset)  
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  subset_obj <- visualizar_resultados(
    subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
    umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
    colores = set_colores, jit_w = 0,
    heatmap = TRUE, heatmap_labels = TRUE, heatmap_fontsize_row = 2,
    volcano_plot = FALSE,
    exp_plot = FALSE, dibujar_violin = FALSE, exp_plot_emparejado = FALSE,
    height = 3600
  )
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}

# EXPRESSION PLOT
for (nombre_subset in c(subsets_para_resultados)) {
  print(nombre_subset)  
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  subset_obj <- visualizar_resultados(
    subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
    umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
    colores = set_colores, jit_w = 0,  n_genes = 10, n_pagina = 1,
    heatmap = FALSE, volcano_plot = FALSE,
    exp_plot = TRUE, dibujar_violin = FALSE, exp_plot_emparejado = FALSE
  )
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}





# Extraer tabla de expresion y guardarla en hoja de calculo
for (nombre_subset in c(subsets_para_resultados)) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funcion
  subset_obj <- guardar_tabla_expresion(subset_obj)
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}





#-----------------------------------------------------------------
#-----------------------------------------------------------------



#////////////////////////////////////////////////////////////
#---------- ANALISIS DE ENRIQUECIMIENTO FUNCIONAL -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# Aplicar a cada subset
for (nombre_subset in subsets_para_resultados) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funciones
  run_enrichment(
    subset_obj = subset_obj,
    label      = nombre_subset,
    OrgDb      = org.Mm.eg.db,
    kegg_org   = "mmu"
  )
}
