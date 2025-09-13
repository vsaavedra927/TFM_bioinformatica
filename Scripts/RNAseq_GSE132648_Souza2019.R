# Resetear
rm(list = ls())


# Establecer ruta de la carpeta
path_carpetaGithub <- "C:/Users/victo/OneDrive/Escritorio/Github"

# Establecer valores de interes
umbral_fdr <- 0.1
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
gse_nr <- "GSE132648"
gse_description <- "GSE132648_Souza2019"
especie <- "Homo sapiens"

# Marcamos la ruta y guardamos nombres de archivo en una lista
path_dataset <- paste0(path_carpetaGithub, "/Datasets/", gse_description, "/")
file_names <- list.files(path = path_dataset, pattern = "\\.txt$")
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
exp_df_all <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

# Convertimos la primera columna en rownames pero asegurando que sean únicos
exp_df_all <- exp_df_all[!duplicated(exp_df_all$X), ]
row.names(exp_df_all) <- exp_df_all$X

# Eliminamos la primera columna (ya está como rownames)
exp_df_all <- exp_df_all[, -1]

# Exploramos rapidamente las primeras filas y la estructura del dataframe de expresion
head(exp_df_all)
str(exp_df_all)




# CREAMOS DF DE GENES CON CODIGOS ENSEMBL (para Gene Ontology)
# (En Souza2019 nos dan simbolos y tenemos que buscar los ids)

# Conectar al dataset humano en Ensembl
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Obtener ENSEMBL IDs a partir de los simbolos de genes
id_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = "external_gene_name",
                         values = row.names(exp_df_all),
                         mart = mart) %>%
            transmute(
              "simbolo gen" = hgnc_symbol,
              "id gen" = ensembl_gene_id
            )




# CREAMOS DATAFRAME DE METADATOS

# Acceso a metadatos usando getGEO
gse <- getGEO(gse_nr, GSEMatrix = TRUE)
names(gse)
names(pData(gse[[1]]))
metadata_df_all <- pData(phenoData(gse[[1]]))

# Extraemos info de tejido y tratamiento en columnas nuevas
metadata_df_all$individuo <- sub("^([0-9]+).*", "\\1", metadata_df_all$title)
metadata_df_all["GSM3884743", "individuo"] <- 3  # Corregimos dato erroneo
metadata_df_all$tratamiento <- metadata_df_all[["sample info:ch1"]]
metadata_df_all$tejido <- metadata_df_all[["tissue:ch1"]]
metadata_df_all <- metadata_df_all %>%
  mutate(
    tratamiento = case_when(
      tratamiento == "Placebo group" ~ "Control",
      tratamiento == "Mineral Oil supplement group" ~ "Aceite",
      TRUE ~ tratamiento   # Deja otros valores como estan
  ),
  tejido = case_when(
    tejido == "Blood samples" ~ "Sangre",
    TRUE ~ tejido
  )
)
# Extraemos los codigos de muestra equivalentes al dataframe de expresion
metadata_df_all <- metadata_df_all %>%
mutate(
  muestra = case_when(
    tratamiento == "Control" ~ paste0("p", title),
    tratamiento == "Aceite"  ~ paste0("o", title),
    TRUE ~ title   # por si hay otros tratamientos
  ),
  muestra = str_replace_all(muestra, "-", "")  # eliminar guiones
)
# Corregimos un codigo de muestra erroneo
metadata_df_all["GSM3884743", "muestra"] <- "o3B24h"


# AJUSTE DE DATOS Y METADATOS

# Nos aseguramos de que los nombres de muestras coinciden en datos y metadatos
print(paste0("Coinciden los nombres de todas las muestras: ", all(metadata_df_all$muestra == colnames(exp_df_all))))
metadata_df_all$muestra == colnames(exp_df_all)
# Hacer columnas de datos igual a filas de metadatos (formato GSM...)
colnames(exp_df_all) <- rownames(metadata_df_all)
print(paste0("Coinciden los nombres de todas las muestras: ", all(rownames(metadata_df_all) == colnames(exp_df_all))))


# CONVERTIR COLUMNAS DE INTERES A FACTOR
# Convertimos variables de interes a factor y marcamos nivel de referencia
metadata_df_all$individuo <- factor(metadata_df_all$individuo)
metadata_df_all$tratamiento <- relevel(factor(metadata_df_all$tratamiento), ref = "Control")
levels(metadata_df_all$tratamiento)


# CONVERTIMOS DATAFRAME A MATRIZ
count_matrix_all <- as.matrix(exp_df_all)


# GENERAR SUBSETS DE DATOS EN FORMATO LISTA
# Generamos lista de propiedades del set completo
set_all <- list(nombre_subset = "all",
                gse_nr = gse_nr,
                especie = especie,
                tejido = "Sangre",
                tratamiento = c("Control", "Aceite"),
                muestras = rownames(metadata_df_all),
                metadata_df = metadata_df_all,
                count_matrix = count_matrix_all
                )

# (Elipsis de codigo para asignar valores de otros subsets)


# COLORES
colores_tratamiento <- setNames(c("#999999", "#E6B800"), unique(metadata_df_all$tratamiento))
# Creamos lista con las paletas de colores
set_colores = list(
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
set_all <- crear_DESeqDataSet(set_all, "individuo + tratamiento")
set_all <-  normalizar_conteos(set_all)
# Obtener matriz de correlacion; graficar heatmap de correlacion y PCA
set_all <- estudiar_correlacion_datos(set_all, colores = set_colores)
# PCA 3D: comprobar el tercer componente principal, en un plot tridimensional
pca_plot_3d(set_all$vsd, set_all$metadata_df, "tratamiento", colores = set_colores$tratamiento, titulo = "PCA 3D - Sangre (Homo sapiens)")





#//////////////////////////////////////////////////////////////
#---- ANALISIS DE EXPRESION DIFERENCIAL (DEA)  CON DESEQ2 -----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# MODELAR EXPRESION DIFERENCIAL, GRAFICAR MEDIA VS VARIANZA
set_all <- aplicar_deseq(set_all)
plot_varvsmean(set_all) 


# CALCULO DE RESULTADOS POR DEFECTO (SIN REDUCCION DE LOG2 FOLD CHANGES)
# Mostrar posible coeficientes de cada subset
resultsNames(set_all$dds)

# Extraer resultados de cada subset, con el coeficiente de interes
# sangre
set_all <- extraer_resultados_log2fcshrink(
                                set_all, 
                                coeficiente = "tratamiento_Aceite_vs_Control", 
                                control_group = "Control",
                                metadatos_global = metadata_df_all,
                                reducir_lfc = FALSE, 
                                umbral_fdr = umbral_fdr,
                                umbral_lfc = umbral_lfc)



#//////////////////////////////////////////////////
#---------- VISUALIZACION DE RESULTADOS -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# Indicar las columnas de metadatos con interes para el expression plot (si se duda, incluir mas columnas)
col_metadata_para_exp_plot <- c("title","geo_accession","organism_ch1",
  "description.1","data_processing.4","sample info:ch1",
  "tissue:ch1","tratamiento","tejido","muestra","individuo")

# Crear plots de visualizacion de los subsets de interes
for (nombre_subset in c("set_all")) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  for (n in 1:1){
    subset_obj <- visualizar_resultados(
                        subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
                        umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
                        colores = set_colores, jit_w = 0.05, n_genes = 15, n_pagina = n,
                        heatmap = TRUE, volcano_plot = TRUE,
                        exp_plot = TRUE, dibujar_violin = TRUE, exp_plot_emparejado = FALSE
                      )
  }
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}


# Heatmap de menos altura
for (nombre_subset in c("set_all")) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  subset_obj <- visualizar_resultados(
    subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
    umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
    colores = set_colores, jit_w = 0.05, n_genes = 15, n_pagina = n,
    heatmap = TRUE, height = 1000, volcano_plot = FALSE,
    exp_plot = FALSE, dibujar_violin = FALSE, exp_plot_emparejado = FALSE
  )
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}


# Solo EXP. PLOT EMPAREJADO para 15 genes
for (nombre_subset in c("set_all")) {
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  subset_obj <- visualizar_resultados(
    subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
    umbral_fdr = umbral_fdr, umbral_lfc = umbral_lfc,
    colores = set_colores, jit_w = 0.05, n_genes = 15,
    heatmap = FALSE, volcano_plot = FALSE,
    exp_plot = FALSE, dibujar_violin = FALSE, exp_plot_emparejado = TRUE
  )
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}




# Extraer tabla de expresion y guardarla en hoja de calculo
for (nombre_subset in c("set_all")) {
  # Recuperar el subset desde el entorno global
  subset_obj <- get(nombre_subset, envir = .GlobalEnv)
  # Aplicar funcion
  subset_obj <- guardar_tabla_expresion(subset_obj)
  # Guardar de nuevo el subset actualizado en el entorno global
  assign(nombre_subset, subset_obj, envir = .GlobalEnv)
}






## -------------------------------------------
## -------------------------------------------

#////////////////////////////////////////////////////////////
#---------- ANALISIS DE ENRIQUECIMIENTO FUNCIONAL -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


run_enrichment(
  set_all, 
  label = "Souza_aceitevsplacebo", 
  OrgDb = org.Hs.eg.db, 
  kegg_org = "hsa")
