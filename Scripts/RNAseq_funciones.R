# --- RNAseq_funciones.R ---

## ------------------------------------------------------------------------
## ////////////////////////////////////////////////////////////////////////
##                   PREPARACION DE ENTORNO Y DATOS
## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## ------------------------------------------------------------------------


#///////////////////
#---- LIBRERIAS ----
#\\\\\\\\\\\\\\\\\\\

.pkgs_bioc <- c(
  "DESeq2",              # DEA
  "SummarizedExperiment",# assay(), vst
  "GEOquery",            # metadatos desde GEO
  "org.Mm.eg.db",        # anotación raton
  "org.Hs.eg.db",        # anotación humano
  "clusterProfiler",     # ORA/GSEA, bitr()
  "enrichplot",          # dotplot, ridgeplot, gseaplot2
  "ReactomePA"           # Reactome: enrichPathway, gsePathway
)

.pkgs_cran <- c(
  "pheatmap",    # heatmaps
  "ggplot2",     # graficos varios
  "dplyr",       # pipes y manipulacion de datos
  "tidyr",       # pivot_longer()
  "tibble",      # rownames_to_column()
  "readr",       # exportar csv
  "stringr",     # str_wrap() en ejes
  "ggrepel",     # etiquetas en volcano
  "patchwork",   # componer graficos
  "plotly",      # PCA 3D
  "scales",      # pretty_breaks() y notacion cientifica
  "ggraph",      # emapplot
  "ggridges"     # ridgeplot
)

load_pkgs <- function(pkgs_bioc = .pkgs_bioc, pkgs_cran = .pkgs_cran){
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
  for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  invisible(lapply(c(pkgs_bioc, pkgs_cran), library, character.only = TRUE))
}



#/////////////////////////////////////
#---- GESTION SUBSETS Y METADATOS ----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# ASIGNAR Nº GSE Y ESPECIE A UN SUBSET
asignar_gse_y_especie <- function(subset_obj, gse_nr, especie) {
  subset_obj$gse_nr <- gse_nr
  subset_obj$especie <- especie
  return(subset_obj)
}


# ASIGNAR MUESTRAS A UN SUBSET
# Funcion para asignar muestras a cada subset relativas a sus propiedades basicas
asignar_muestras <- function(subset_obj, metadatos_global) {
  # Filtrar metadata_df_all y obtener geo_accession
  subset_obj$muestras <- metadatos_global %>%
    dplyr::filter(tejido == subset_obj$tejido, tratamiento %in% subset_obj$tratamiento) %>%
    dplyr::pull(geo_accession)
  return(subset_obj)
}


# ASIGNAR METADATOS A UN SUBSET
# Funcion para asignar metadatos correspondientes a cada subset
asignar_metadatos <- function(subset_obj, metadatos_global){
  # Generar subset de 'metadatos'
  subset_obj$metadata_df <- metadatos_global[metadatos_global$geo_accession %in% subset_obj$muestras, ]
  message(paste0("Objeto creado:   ", subset_obj$nombre_subset, "$metadata_df"))
  return(subset_obj)
}


# ASIGNAR MATRIZ DE CONTEOS A UN SUBSET
# Funcion para asignar matriz de conteos correspondientes a cada subset
asignar_matriz_conteos <- function(subset_obj, count_matrix_global){
  # Generar subset de 'metadatos'
  subset_obj$count_matrix <- count_matrix_global[,subset_obj$muestras]
  message(paste0("Objeto creado:   ", subset_obj$nombre_subset, "$count_matrix"))
  return(subset_obj)
}


# CREAR COPIAS DE UN SUBSET CON SUFIJO
copiar_subset <- function(subset_obj, sufijo){
  # Modificar el campo nombre_subset
  subset_obj$nombre_subset <- paste0(subset_obj$nombre_subset, sufijo)
  message(paste0("Creado nuevo subset: ", subset_obj$nombre_subset))
  return(subset_obj)
}


## ------------------------------------------------------------------------
## ////////////////////////////////////////////////////////////////////////
##                ANALISIS DE EXPRESION DIFERENCIAL
## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## ------------------------------------------------------------------------


#///////////////////////////////////
#--------- ANALISIS DESEQ2 ---------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# OBJETO DESEQDATASET
# Funcion para crear deseqdataset correspondiente a un subset
crear_DESeqDataSet <- function(subset_obj, variable_diseno){
  # Convertimos el nombre de la variable a formula
  formula_diseno <- as.formula(paste0("~ ", variable_diseno))
  # Creamos el objeto DESeqDataSet: dds
  subset_obj$dds <- DESeqDataSetFromMatrix(countData = subset_obj$count_matrix,
                                           colData = subset_obj$metadata_df,
                                           design = formula_diseno)
  message(paste0("Objeto creado:   ", subset_obj$nombre_subset, "$dds"))
  return(subset_obj)
}            



# NORMALIZACION DE LOS CONTEOS
# Funcion para normalizar los conteos de un subset
normalizar_conteos <- function(subset_obj){
  # Determinar los size factors a usar en la normalizacion, en el DESeqDataSet o dds
  subset_obj$dds <- estimateSizeFactors(subset_obj$dds)
  print(paste0("Estimados los factores de tamaño de ", subset_obj$nombre_subset))
  # Extraer el dataframe de conteos normalizados
  subset_obj$count_matrix_norm <- counts(subset_obj$dds, normalized=TRUE)
  message(paste0("Objeto creado:   ", subset_obj$nombre_subset, "$count_matrix_norm"))
  return(subset_obj)
}



#///////////////////////////////////
#---- EXPLORACION DE LOS DATOS -----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# EXPLORAR CORRELACION DE MUESTRAS
# 1. CALCULAR MATRIZ DE CORRELACION
# 2. CLUSTERIZACION NO SUPERVISADA -> HEATMAP DE CORRELACION (opcional)
# 3. CLUSTERIZACION NO SUPERVISADA -> PRINCIPAL COMPONENT ANALYSIS (PCA) EN 2D (opcional)
# La lectura de la correlacion nos permite:
# · Hacer un analisis de calidad de las muestras (se agrupan muestras similares entre si como esperamos? hay outliers?)
# ·Encontrar fuentes adicionales de variacion entre las variables (clusterizacion mas fuerte en unas variables que en otras: cepa, sexo, tratamiento...)
#    (Las variables encontradas en 2. se usaran en la design formula en los pasos siguientes)

estudiar_correlacion_datos <- function(subset_obj, 
                                       heatmap = TRUE, pca = TRUE, 
                                       colores = NULL, file_name = NULL, 
                                       width = 1200, height = 1000, res = 150){
  
  # 1) CALCULAR MATRIZ DE CORRELACION DE MUESTRAS, si no existe aun
  # Aplicar transformacion de varianza estabilizada
  # (aproxima conteos RNA-seq con distribucion binomial negativa a una distribucion normal)
  if (is.null(subset_obj$vsd_cor)){
    print(paste0("Aplicando transformacion de varianza estabilizada (vst) de", subset_obj$nombre_subset))
    subset_obj$vsd <- DESeq2::vst(subset_obj$dds, blind = TRUE)
    # Extraer la matriz vst del objeto (genes x muestras)
    print(paste0("Extrayendo matriz vst de ", subset_obj$nombre_subset))
    subset_obj$vsd_mat <- SummarizedExperiment::assay(subset_obj$vsd)
    # Calcular matriz de correlacion entre muestras (correlacion de Pearson)
    print(paste0("Calculando matriz de correlacion de muestras de ", subset_obj$nombre_subset))
    subset_obj$vsd_cor <- stats::cor(subset_obj$vsd_mat)
  }
  
  # 2) HEATMAP DE CORRELACION, opcional
  if (heatmap){
    # Indicar título y ruta de archivo
    etiqueta <- paste0("Mapa de calor de correlación de muestras\n", subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", paste(subset_obj$tratamiento, collapse = ", "))
    file_name <- paste0("correlation_heatmap_", subset_obj$nombre_subset , ".png")
    # Generar plot y guardarlo como archivo png
    # annotation: que factor de metadatos incluir como barras de anotacion
    pheatmap(subset_obj$vsd_cor, 
             annotation = dplyr::select(subset_obj$metadata_df, tratamiento), 
             annotation_colors = colores,
             main = etiqueta,
             filename = file.path(path_results, file_name),
             width = width/res, height = height/res   # pheatmap trabaja en pulgadas; pulgadas = pixeles / resolucion
    )
    message(paste0("Guardado heatmap de correlacion: ", file_name))
  }
  
  # 3) PRINCIPAL COMPONENT ANALYSIS (PCA), opcional
  if(pca){
    # Generar plot
    titulo <- "Análisis de componentes principales (PCA)"
    subtitulo <- paste0(subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", paste(subset_obj$tratamiento, collapse = ", "))
    pca_plot <- plotPCA(subset_obj$vsd, intgroup = "tratamiento") +
      ggtitle(titulo, subtitle = subtitulo) +
      labs(color = "tratamiento") +
      scale_color_manual(values = colores$tratamiento)
    # Guardar en archivo
    file_name <- paste0("pca_plot_", subset_obj$nombre_subset , ".png")
    ggplot2::ggsave(
      filename = file.path(path_results, file_name),
      plot = pca_plot,
      width = 6,
      height = 5,
      dpi = 300
    )
    message(paste0("Guardado PCA: ", file_name))
  }
  return(subset_obj)
}



# GRAFICAR PCA 3D
# (Requiere haber obtenido objeto vsd de la funcion anterior)
pca_plot_3d <- function(vsd, metadata, grupo, colores = NULL, titulo = "PCA 3D") {
  # Realizar PCA con prcomp
  pca <- prcomp(t(SummarizedExperiment::assay(vsd)))
  # Porcentaje de varianza explicada
  var_exp <- round(100 * summary(pca)$importance[2, 1:3], 1)
  # Extraer las tres primeras PCs
  pca_df <- as.data.frame(pca$x[, 1:3])
  pca_df[[grupo]] <- metadata[[grupo]]
  # Crear gráfico interactivo con plotly
  plot_ly(pca_df, 
          x = ~PC1, y = ~PC2, z = ~PC3, 
          color = ~get(grupo),
          colors = colores,
          type = "scatter3d", 
          mode = "markers") %>%
    layout(title = titulo,
           scene = list(
             xaxis = list(title = paste0("PC1 (", var_exp[1], "%)")),
             yaxis = list(title = paste0("PC2 (", var_exp[2], "%)")),
             zaxis = list(title = paste0("PC3 (", var_exp[3], "%)"))
           ))
}




#//////////////////////////////////////////////////////////////
#---- ANALISIS DE EXPRESION DIFERENCIAL (DEA)  CON DESEQ2 -----
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# MODELAR EXPRESION DIFERENCIAL (segun distribucion binomial negativa) PARA CADA GEN
# Aplicamos DESeq a los DESeqDataSet
# DESeq ajusta un modelo lineal generalizado de tipo binomial negativa usando la formula de diseño
# Pasos internos de DESeq:
# 1. Estimacion de factores de tamaño, estimateSizeFactors
# 2. Estimacion de dispersion, estimateDispersions
# 3. Ajuste a Binomial Negativa y test estadistico de Wald, nbinomWaldTest
aplicar_deseq <- function(subset_obj){
  # Aplicar DESeq al subset
  subset_obj$dds <- DESeq(subset_obj$dds)
  message(paste0("Aplicado DESeq a:   ", subset_obj$nombre_subset))
  return(subset_obj)
}



# PLOT DE VARIANZA VS MEDIA
# Permite comprobar que la distribucion se comporta como binomial negativa
plot_varvsmean <- function(subset_obj){
  # Calcular media y varianzas de los conteos para cada gen
  mean_counts <- apply(subset_obj$count_matrix, 1, mean)
  variance_counts <- apply(subset_obj$count_matrix, 1, var)
  # Crear dataframe con media y varianza de cada gen
  subset_obj$meanvar_df <- data.frame(
    mean = mean_counts, 
    var = variance_counts, 
    row.names = rownames(subset_obj$count_matrix)
  )
  # Plotear
  titulo <- "Media vs varianza"
  subtitulo <- paste0(subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", paste(subset_obj$tratamiento, collapse = ", "))
  meanvar_plot <- ggplot(subset_obj$meanvar_df) +
    geom_point(aes(x=mean, y=var)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Línea Poisson
    scale_y_log10() +
    scale_x_log10() +
    xlab("Media por gen") +
    ylab("Varianza por gen") +
    ggtitle(titulo, subtitle = subtitulo)
  message(paste0("Creado plot varianza vs media:   ", subset_obj$nombre_subset, "$meanvar_plot"))
  # Guardar
  file_name = paste0("meanvar_plot_", subset_obj$nombre_subset , ".png")
  ggplot2::ggsave(
    filename = file.path(path_results, file_name),
    plot = meanvar_plot,
    width = 6,
    height = 5,
    dpi = 300
  )
  print(paste0("Guardado plot Varianza vs Media:   ", file_name))
  # Generar plot de estimaciones de dispersion
  etiqueta <- paste0("Estimaciones de dispersión\n", subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", paste(subset_obj$tratamiento, collapse = ", "))
  file_name = paste0("estimateddispersion_plot_", subset_obj$nombre_subset , ".png")
  png(filename = file.path(path_results, file_name), width = 6, height = 5, units = "in", res = 300)
  plotDispEsts(subset_obj$dds, main = etiqueta)
  message(paste0("Guardado plot de estimacion de dispersiones:   ", file_name))
  dev.off()
}




# CALCULO DE RESULTADOS POR DEFECTO (A) O CON REDUCCION DE LOG2 FOLD CHANGES (B)
# results() extrae del objeto dds un resumen estadistico para todas las comparaciones definidas por el diseño
extraer_resultados_log2fcshrink <- function(subset_obj, 
                                            coeficiente, 
                                            metadatos_global,
                                            control_group = "Control",
                                            reducir_lfc = FALSE, 
                                            umbral_fdr = 0.05,
                                            umbral_lfc = 0.585){
  message(paste0("Procediendo a extraer resultados de:   ", subset_obj$nombre_subset, 
                 "\nReducción de log2FC: ", reducir_lfc,
                 "\nUmbral de FDR para genes significativos: ", umbral_fdr,
                 "\nUmbral de log2FC para genes significativos: ", umbral_lfc))
  # Comprobar que el coeficiente es correcto
  if (!(coeficiente %in% resultsNames(subset_obj$dds))) {stop("No se ha encontrado el coeficiente")}
  
  # Indicar constantes para ambos tipos (iguales para resultados con o sin reduccion de LFC)
  deteccion_outliers<- TRUE
  filtrado_independiente <- TRUE
  
  # 1(A). Extraer resultados de exp.dif. por defecto, filtrando genes por padj < 0.05
  subset_obj$results <- results(subset_obj$dds, 
                                name = coeficiente,
                                alpha = umbral_fdr,
                                independentFiltering = filtrado_independiente,
                                cooksCutoff = deteccion_outliers
                                )
  message(paste0("Extraidos resultados por defecto en:   ", subset_obj$nombre_subset, "$results\n"))
  if (!reducir_lfc) {
    # EXPLORACION DE RESULTADOS sin reduccion de LFC
    # Mostrar tabla de resultados para los 10 primeros genes, y descripcion de las columnas
    print(as.data.frame(mcols(subset_obj$results)))
    #message(paste0("Extraidos resultados con reducción de log2FC en:   ", subset_obj$nombre_subset, "$results_lfcshrink")
    cat("Genes sin ordenar ni filtrar. Se muestran solo los 10 primeros genes)\n-\n")
    print(head(subset_obj$results, n=10))
    # Extraer y ordenar los genes significativos (p-value ajustado < umbral_fdr (por defecto 0.05))
    subset_obj$results_df <- data.frame(subset_obj$results)
    rownames(subset_obj$results_df) <- rownames(subset_obj$results)
  }
  
  # 1(B). Extraer resultados de exp.dif. con log2FC shrinkage, para conseguir log2FC mas precisos
  # Esto no cambia los p-values ni los p-adj, por tanto no deberia modificar que genes son significativos
  # Version con "apeglm", menos sesgado
  if (reducir_lfc) {
    subset_obj$results_lfcshrink <- lfcShrink(subset_obj$dds, 
                                              coef = coeficiente,
                                              res = subset_obj$results,
                                              type = "apeglm")
    subset_obj$results_lfcshrink
    message(paste0("Extraidos resultados con reducción de log2FC en:   ", subset_obj$nombre_subset, "$results_lfcshrink",
                   "\n-"))
    # EXPLORACION DE RESULTADOS tras reduccion de LFC
    # Mostrar tabla de resultados para los 10 primeros genes, y descripcion de las columnas
    print(as.data.frame(mcols(subset_obj$results_lfcshrink)))
    #message(paste0("Extraidos resultados con reducción de log2FC en:   ", subset_obj$nombre_subset, "$results_lfcshrink")
    cat("Genes sin ordenar ni filtrar. Se muestran solo los 10 primeros genes)\n-\n")
    print(head(subset_obj$results_lfcshrink, n=10))
    # Extraer y ordenar los genes significativos (p-value ajustado < umbral_fdr (por defecto 0.05))
    subset_obj$results_df <- data.frame(subset_obj$results_lfcshrink)
    rownames(subset_obj$results_df) <- rownames(subset_obj$results_lfcshrink)
  }
  
  # 2. ESTRELLAS DE SIGNIFICACION
  # Pasar los genes a columnas (columna a primera posicion)
  subset_obj$results_df$gene <- rownames(subset_obj$results_df)
  subset_obj$results_df <- subset_obj$results_df %>%
    dplyr::relocate(gene) 
  # Ordenar los genes por padj
  subset_obj$results_df <- subset_obj$results_df %>%
    arrange(padj)
  # Añadir estrellas de significancia
  subset_obj$results_df <- subset_obj$results_df %>%
    mutate(
      sig_stars = case_when(
        is.na(padj)        ~ "ns",
        padj < 1e-4        ~ "****",
        padj < 1e-3        ~ "***",
        padj < 1e-2        ~ "**",
        padj < umbral_fdr  ~ "*",
        TRUE            ~ "ns"   # Otros valores
      ))
  # Incluir el valor del máximo del grupo de tratamiento (no control): asi calcularemos posicion de estrellas
  muestras_tratamiento_global <- rownames(metadatos_global)[metadatos_global$tratamiento %in% setdiff(subset_obj$tratamiento, control_group)]
  muestras_tratamiento_subset <- subset_obj$muestras[subset_obj$muestras %in% muestras_tratamiento_global]
  subset_obj$results_df <- subset_obj$results_df %>%
    mutate(
      y_max = apply(
        subset_obj$count_matrix_norm[
          rownames(subset_obj$results_df),  # Usamos el orden de filas del df de resultados
          muestras_tratamiento_subset,   # Solo nos interesan las muestras de tratamiento (no control)
          drop = FALSE], 
        1, max, na.rm = TRUE))
  # Extraer genes significativos (p-value ajustado < umbral_fdr (por defecto 0.05), |log2FC| > umbral_lfc (por defecto 0.585))
  subset_obj$results_df_sig_neg <- subset_obj$results_df %>%
                                      filter(padj < umbral_fdr) %>%
                                      filter(log2FoldChange < (0 - umbral_lfc))
  subset_obj$results_df_sig_pos <- subset_obj$results_df %>%
                                      filter(padj < umbral_fdr) %>%
                                      filter(log2FoldChange > umbral_lfc)
  subset_obj$results_df_sig <- bind_rows(subset_obj$results_df_sig_pos, subset_obj$results_df_sig_neg) %>%
                                      arrange(padj)
  # 3. Mostrar resultados
  cat("\n--- GENES SIGNIFICATIVOS ---\n")
  n_sig_pos <- nrow(subset_obj$results_df_sig_pos)
  n_sig_neg <- nrow(subset_obj$results_df_sig_neg)
  cat(n_sig_pos, " genes sobre-expresados significativamente (log2FC > ", umbral_lfc, ", padj < ", umbral_fdr, "):  ", paste(rownames(subset_obj$results_df_sig_pos), collapse = ", "), "\n",
      n_sig_neg, " genes infra-expresados significativamente (log2FC < ", umbral_lfc, ", padj < ", umbral_fdr, "):  ", paste(rownames(subset_obj$results_df_sig_neg), collapse = ", "), "\n", sep="")
  cat(" \n")
  cat("Total de genes significativos: ", nrow(subset_obj$results_df_sig), "\n-")
  message(paste0("\nGenes significativos almacenados en:   ", subset_obj$nombre_subset, "$results_df_sig"))
  summary(subset_obj$results, alpha = umbral_fdr)
  
  # 4. Si hay genes significativos, mostrar un resumen y calcular conteos normalizados para ellos
  if (nrow(subset_obj$results_df_sig) > 0){
    cat("Mostrando genes más significativos (menor valor p ajustado):\n-\n")
    print(head(subset_obj$results_df_sig, 10))
    cat("...\n")
    # Generar subset de los conteos normalizados con solo los genes significativos, util para el expression plot
    subset_obj$count_matrix_norm_sig <- subset_obj$count_matrix_norm[rownames(subset_obj$results_df_sig), ]
    message(paste0("Generada matriz de conteos normalizados para los genes significativos en:   ", subset_obj$nombre_subset, "$count_matrix_norm_sig"))
    # Mostrar primeras filas de count_matrix_norm_sig
    head(subset_obj$count_matrix_norm_sig, 5)
    cat("...\n")
  } else {
    message(paste0("No existen genes significativos en   ", subset_obj$nombre_subset, "$results_df_sig"))
    cat(" \n")
  }
  return(subset_obj)
}




#//////////////////////////////////////////////////
#---------- VISUALIZACION DE RESULTADOS -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# GRAFICOS DE VISUALIZACION DE RESULTADOS
# A) HEATMAP DE GENES SIGNIFICATIVOS
# B) VOLCANO PLOT (con etiquetas de genes significativos) DE TODOS LOS GENES
# C) JITTER PLOT ("expression plot") DE CADA GEN SIGNIFICATIVO, con opcion de VIOLIN PLOT
# D) JITTER PLOT EMPAREJADO DE CADA GEN SIGNIFICATIVO
# Nota 1: Para mayor control de la visualizacion, cada grafico es opcional.
# Nota 2: si se activa un jitter plot, se calcula un df adicional en formato largo

visualizar_resultados <- function(subset_obj, metadatos_global, col_metadata, control_group = "Control", 
                                  umbral_lfc = 0.585, umbral_fdr = 0.05,
                                  colores = NULL, jit_w = 0, n_genes = 10, n_pagina = 1,
                                  heatmap = TRUE, heatmap_labels = TRUE, heatmap_fontsize_row = 6,
                                  volcano_plot = TRUE, 
                                  exp_plot = TRUE, dibujar_violin = TRUE, exp_plot_emparejado = TRUE, 
                                  width = 1200, height = 1200, res = 150){
  
  # A) EXPRESSION HEATMAP: ploteamos los conteos normalizados de los genes significativos
  if (heatmap){
    # Saltar este plot si no hay genes significativos para representar
    if (nrow(subset_obj$results_df_sig) == 0) {
      message("El subset ", subset_obj$nombre_subset, " no tiene genes significativos para representar en el heatmap.")
    } else {
      # Indicar titulo y ruta de archivo
      etiqueta = paste0("Mapa de calor de expresión de genes significativos\n", 
                        subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", 
                        paste(subset_obj$tratamiento, collapse = ", "), 
                        "\nFDR < ", umbral_fdr, " ; |log2FC| > ", umbral_lfc)
      file_name = paste0("expression_heatmap_", subset_obj$nombre_subset , ".png")
      # Generar heatmap y guardarlo como archivo png
      pheatmap(subset_obj$count_matrix_norm_sig,
               color = colorRampPalette(c("blue", "white", "red"))(100),  #   gradiente de colores: azul (bajo) → blanco (0) → rojo (alto)
               breaks = seq(-2, 2, length.out = 101),   # truncar valores extremos de color, más resolución en torno a 0
               cluster_rows = TRUE,
               show_rownames = heatmap_labels, 
               fontsize_row = heatmap_fontsize_row,
               show_colnames = TRUE,
               annotation = dplyr::select(subset_obj$metadata_df, tratamiento),  # barras de anotacion, que factor de metadatos mostrar
               annotation_colors = colores,
               scale = "row",
               display_numbers = FALSE,
               main = etiqueta,
               filename = file.path(path_results, file_name),
               width = width/res, height = height/res
               )
      message(paste0("Guardado heatmap de expresion de genes significativos: ", file_name))
      }
    }
  
  # B) VOLCANO PLOT
  if (volcano_plot){
    # Marcar genes significativos, con p-values ajustados < umbral_fdr (por defecto 0.05), en nueva columna "threshold"
    subset_obj$results_df <- subset_obj$results_df %>%
      mutate(threshold = padj < umbral_fdr) 
    # Título del plot
    titulo <- "Volcano plot"
    subtitulo <- paste0(subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", 
                        paste(subset_obj$tratamiento, collapse = ", "), 
                        "\nFDR < ", umbral_fdr, " ; |log2FC| > ", umbral_lfc)
    
    # Composicion en patchwork del volcano plot y la leyenda.
    # 1. Volcano plot
    subset_obj$volcano_plot <- ggplot(subset_obj$results_df, aes(log2FoldChange,-log10(pvalue))) + 
      # Grafico de puntos (cambiado a gris mas claro); nombres de ejes x, y
      geom_point(colour = "gray", size = 3) +
      scale_y_continuous(name = "-log10 p-value") +
      scale_x_continuous(name = "log2 fold change") +
      ggtitle(titulo, subtitle = subtitulo) +
      # Colorear puntos significativos
      # · ejecutar solo si hay genes significativos
      # · Genes sobreexpresados en rojo, infraexpresados en azul
      {if (!is.null(subset_obj$results_df_sig_pos) && nrow(subset_obj$results_df_sig_pos) > 0) 
        geom_point(data = subset_obj$results_df_sig_pos, 
                   aes(log2FoldChange,-log10(pvalue)), color = "red3", size = 4) } + 
      {if (!is.null(subset_obj$results_df_sig_neg) && nrow(subset_obj$results_df_sig_neg) > 0) 
        geom_point(data = subset_obj$results_df_sig_neg, 
                   aes(log2FoldChange,-log10(pvalue)), color = "royalblue3", size = 4) } +
      # Personalizacion: sin leyenda; tamaño y color de ejes; sin cuadricula (grid) de fondo; margenes amplios
      theme(legend.position = "none",
            axis.title = element_text(size = 25),
            axis.text.x  =  element_text(size = 20, hjust = 1, colour = "black"),
            axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"),
            plot.title = element_text(size = 40, hjust = 0.5),   
            plot.subtitle = element_text(size = 28, hjust = 0.5),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.margin=unit(c(0, 0, 1, 3), "cm"))
    
    # 2. Recuadro para leyenda
    leyenda <- ggplot() +
      xlim(0,1) + ylim(0,1) +
      annotate("rect", xmin=0, xmax=1, ymin=0.65, ymax=0.35,
               fill="white", colour="black") +
      annotate("point", x=.07, y=.6, colour="royalblue3", size=4) +
      annotate("text",  x=.17, y=.6, label="Infraexpresado", hjust=0, size=8) +
      annotate("point", x=.07, y=.50, colour="grey50", size=3) +
      annotate("text",  x=.17, y=.50, label="No significativo", hjust=0, size=8) +
      annotate("point", x=.07, y=.4, colour="red3", size=4) +
      annotate("text",  x=.17, y=.4, label="Sobreexpresado", hjust=0, size=8) +
      theme_void() +
      theme(plot.margin = margin(0,0,0,0))   # un poco de aire 
    
    # 3. Composicion de patchwork (volcano + leyenda)
    final_plot <- subset_obj$volcano_plot + leyenda + plot_layout(widths = c(1, 0.25))
    
    # 4. Guardar composicion de volcano plot con leyenda
    file_name = paste0("volcano_plot_", subset_obj$nombre_subset , ".png")
    ggplot2::ggsave(
      filename = file.path(path_results, file_name),
      plot = final_plot,
      width = 25,
      height = 12,
      dpi = 300, 
      scale = 0.8
    )
    
    message(paste0("Guardado volcano plot: ", file_name))
    
    # Si existen genes significativos
    if (nrow(subset_obj$results_df_sig) > 0){
      # Guardar otro volcano plot con etiquetas de genes significativos
      subset_obj$volcano_plot_labels <- subset_obj$volcano_plot +
        # Etiquetas de genes sobreexpresados
        {if (!is.null(subset_obj$results_df_sig_pos) && nrow(subset_obj$results_df_sig_pos) > 0) 
          geom_text_repel(
            data = subset_obj$results_df_sig_pos, 
            aes(label = gene), 
            size = 8,
            box.padding   = unit(0.6, "lines"),
            point.padding = unit(0.35, "lines"),
            force = 1.6,
            force_pull = 0.6,
            nudge_y = 0.05,        # pequeño empujon inicial (en -log10 p-value)}
            min.segment.length = 2,
            segment.color = "grey50"
            )} +
        # Etiquetas de genes infraexpresados
        {if (!is.null(subset_obj$results_df_sig_neg) && nrow(subset_obj$results_df_sig_neg) > 0) 
          geom_text_repel(
            data = subset_obj$results_df_sig_neg, 
            aes(label = gene), 
            size = 8,
            box.padding   = unit(0.6, "lines"),
            point.padding = unit(0.35, "lines"),
            force = 1.6,
            force_pull = 0.6,
            nudge_y = 0.05,        # pequeño empujon inicial (en -log10 p-value)}
            min.segment.length = 2,
            segment.color = "grey50"
            )}
      # Composicion de patchwork (volcano con etiquetas + leyenda)
      final_plot <- subset_obj$volcano_plot_labels + leyenda + plot_layout(widths = c(1, 0.25))
      # Guardar composicion
      file_name = paste0("volcano_plot_labels_", subset_obj$nombre_subset , ".png")
      ggplot2::ggsave(
        filename = file.path(path_results, file_name),
        plot = final_plot,
        width = 25,
        height = 12,
        dpi = 300, 
        scale = 0.8
      )
      message(paste0("Guardado volcano plot con etiquetas: ", file_name))
    }
  }
  
  
  # (!) PREPARAR TOP GENES EN FORMATO LARGO para expression plots
  if (exp_plot || exp_plot_emparejado){
    # Saltar esto si no hay genes significativos para representar
    if (nrow(subset_obj$results_df_sig) == 0) {
      message("El subset ", subset_obj$nombre_subset, " no tiene genes significativos")
    } else {
      # Definir el numero de genes a incluir en el plot; máx 20 genes, o el total de genes si hay menos de 20  
      nr_top <- min(n_genes, nrow(subset_obj$count_matrix_norm_sig))
      # El parametro n_pagina permite graficar a lo largo de la lista de genes (en caso de numero alto de genes)
      primer_gen_a_graficar <- (n_pagina - 1) * nr_top + 1
      ultimo_gen_a_graficar <- n_pagina * nr_top
      # Convertir matriz de conteos normalizados de genes significativos a un dataframe de formato largo
      subset_obj$top_genes_longdf <-  data.frame(subset_obj$count_matrix_norm_sig)[primer_gen_a_graficar:ultimo_gen_a_graficar,] %>%
        tibble::rownames_to_column("genes") %>%
        # Formato largo: reducimos todas las columnas (genes, muestra1, muestra2, ... muestraN) a tres:
        # · genes: los nombres de los genes
        # · samplename: aqui se reunen todas las columnas de muestras en una sola
        # · normalized_counts: aqui van todos los valores numericos
        tidyr::pivot_longer(cols = -genes,
                            names_to = "samplename",
                            values_to = "normalized_counts") %>%
        # Unir la metadata para colorear muestras por grupos
        dplyr::inner_join(
          metadatos_global[, col_metadata],
          by = c("samplename" = "geo_accession")
        )
      # Convertimos los genes a factor (necesario para plot de violin)
      subset_obj$top_genes_longdf$genes <- factor(subset_obj$top_genes_longdf$genes, 
                                                  levels = unique(subset_obj$top_genes_longdf$genes))
      
      # Calcular n para cada grupo. Texto "Control n=4; DHAsup n=4; EPAsup n=4" (solo los que existan en este subset)
      n_por_grupo <- table(factor(subset_obj$metadata_df$tratamiento, 
                                  levels = subset_obj$tratamiento))
      
      # Para ESTRELLAS DE SIGNIFICANCIA: data frame con posiciones y etiquetas
      # Indicar nombres de genes a representar y niveles de tratamiento
      genes_top <- levels(subset_obj$top_genes_longdf$genes)
      trat_levels <- levels(droplevels((subset_obj$metadata_df$tratamiento)))
      trat_treated <- intersect(setdiff(trat_levels, control_group), subset_obj$tratamiento)                       
      # Generar df con posiciones y etiquetas del tratamiento
      subset_obj$estrellas_df <- subset_obj$results_df_sig %>%
        dplyr::filter(gene %in% genes_top) %>%    # solo los genes del plot
        dplyr::transmute(
          genes = factor(gene, levels = genes_top),   # alinear con el eje x del plot
          ypos = y_max * 2,            # un poco por encima del máximo
          label = sig_stars
        ) 
      subset_obj$estrellas_df <- dplyr::bind_rows(
        subset_obj$estrellas_df,
        # fila vacia (con "") para el control para que el dodge alinee a derecha/izquierda
        subset_obj$estrellas_df %>% dplyr::mutate(
                                              tratamiento = factor(control_group, levels = trat_levels),
                                              label = "")
      )
    }
  }
  
  
  # 3) JITTER PLOT ("EXPRESSION PLOT")
  if (exp_plot) {
    # Saltar este plot si no hay genes significativos para representar
    if (nrow(subset_obj$results_df_sig) == 0) {
      message("El subset ", subset_obj$nombre_subset, " no tiene genes significativos para representar en el plot de violines")
    } else {
      # Indicar titulo y subtitulo
      titulo = paste0("Niveles de expresión normalizados de los genes\ndiferencialmente expresados más significativos")
      subtitulo = paste0(nr_top, " genes (", primer_gen_a_graficar, " a ", ultimo_gen_a_graficar , ") de un total de ", nrow(subset_obj$count_matrix_norm_sig), " genes significativos\n",
                         subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", 
                         paste(paste0(names(n_por_grupo), " (n=", n_por_grupo, ")"), collapse = ", "), 
                         "\nFDR < ", umbral_fdr, " ; |log2FC| > ", umbral_lfc)
      # Indicar ancho de jitter
      dodge_w <- 0.6
      #jit_w indicado en parámetros
      
      # Generar expression plot
      subset_obj$exp_plot <- ggplot(subset_obj$top_genes_longdf) +
        # VIOLINES de distribución
        {if (dibujar_violin) geom_violin(
          aes(x = genes, y = normalized_counts, fill = tratamiento,
              group = interaction(genes, tratamiento)),
          position = position_dodge(dodge_w), 
          width = .85, alpha = .25, trim = FALSE, linewidth = 0.2,
          show.legend = FALSE
        )} +
        # PUNTOS de datos
        geom_point(
          aes(x = genes, y = normalized_counts, color = tratamiento),   # Valores de eje x, eje y, color
          position = position_jitterdodge(jitter.width = jit_w, dodge.width = dodge_w),   # Efecto jitter o desplazamiento de puntos
          alpha = .7, size = 1
        ) +
        # MEDIAS de cada grupo de datos
        stat_summary(
          aes(x = genes, y = normalized_counts, group = interaction(genes, tratamiento)),
          fun = mean, geom = "point", shape = 95, size = 6,
          position = position_dodge(dodge_w), color = "grey20"
        ) +
        # ESTRELLAS DE SIGNIFICANCIA
        geom_text(
          data = subset_obj$estrellas_df,
          aes(x = genes, y = ypos, label = label,
              group = interaction(genes, tratamiento)),
          position = position_dodge(width = dodge_w),
          inherit.aes = FALSE,
          size = 5,
          color = "grey20"
        ) +
        # EJES x,y
        scale_y_log10() +
        xlab("Genes") + ylab("Conteos normalizados") +
        scale_x_discrete(expand = expansion(add = 0.6)) +
        # COLORES
        scale_fill_manual(name = "tratamiento", values = colores$tratamiento) +  # Colores de violin
        scale_color_manual(name = "tratamiento", values = colores$tratamiento) +  # Colores de puntos
        guides(color = guide_legend(override.aes = list(size = 2))) +   # Guia para los colores de leyenda
        # TITULO y TEMAS de fondo
        ggtitle(titulo, subtitle = subtitulo) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              legend.key = element_rect(fill = "white"),
              legend.background = element_rect(fill = "white", color = "white")
        )
      
      # Guardar plot
      # ancho proporcional (1.5" por gen); evita que quede demasiado estrecho
      n_genes <- length(levels(subset_obj$top_genes_longdf$genes))
      width_per_gene <- 1.2
      legend_w   <- 2.5   
      
      width_in   <- n_genes * width_per_gene + legend_w
      
      genes_a_graficar <- paste0(primer_gen_a_graficar, "a", ultimo_gen_a_graficar)
      file_name = paste0("expression_plot_", subset_obj$nombre_subset , "_", genes_a_graficar, ".png")
      ggplot2::ggsave(
        filename = file.path(path_results, file_name),
        plot = subset_obj$exp_plot,
        width = width_in, height = 5, dpi = 300, scale = 0.8
      )
      message(paste0("Guardado gráfico de violines con la expresion de genes significativos: ", file_name))
    }
  }
  
  # 4) JITTER PLOT EMPAREJADO ("EXPRESSION PLOT EMPAREJADO")
  if (exp_plot_emparejado) {
    # Saltar este plot si no hay genes significativos para representar
    if (nrow(subset_obj$results_df_sig) == 0) {
      message("El subset ", subset_obj$nombre_subset, " no tiene genes significativos para representar en el plot de expresión emparejado")
    } else {
      # Indicar titulo y subtitulo
      titulo = paste0("Niveles de expresión normalizados de los genes\ndiferencialmente expresados más significativos")
      subtitulo = paste0(nr_top, " genes de un total de ", nrow(subset_obj$count_matrix_norm_sig), " genes significativos\n",
                         subset_obj$gse_nr, " (", subset_obj$especie, ") - ", subset_obj$tejido, "\n", 
                         paste(subset_obj$tratamiento, collapse = ", "), 
                         "\nFDR < ", umbral_fdr, " ; |log2FC| > ", umbral_lfc)
      # Indicar ancho de jitter
      dodge_w <- 0.6
      #jit_w establecido como parametro de funcion

      # Generar plot
      subset_obj$exp_plot_emparejado <- ggplot(
        subset_obj$top_genes_longdf, 
        aes(x = tratamiento, y = normalized_counts)) +
        ## LÍNEAS DE EMPAREJAMIENTO: agrupa por gen e individuo
        geom_line(
          aes(group = interaction(genes, individuo)),
          color = "grey70", linewidth = 0.4, alpha = 0.7
          ) +
        # PUNTOS de datos
        geom_point(
          aes(color = tratamiento), 
          size = 2,
          position = position_jitter(width = jit_w, height = 0)
          ) +
        # MEDIAS de cada grupo de datos
        stat_summary(
          aes(color = tratamiento),
          fun = mean, geom = "point", shape = 95, size = 12,
          position = position_dodge(dodge_w), color = "grey20"
          ) +
        # PANELES
        facet_wrap(~genes, scales = "free_y") +   # un panel por gen
        # COLORES
        scale_color_manual(values = colores$tratamiento) +
        # TITULOS y TEMAS de fondo
        ggtitle(titulo, subtitle = subtitulo) +
        theme_bw((base_size = 12)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      # Guardar plot
      n_genes <- length(unique(subset_obj$top_genes_longdf$genes))
      n_rows  <- ceiling(n_genes / 4)   # asumiendo 4 columnas en facet_wrap
      file_name = paste0("expression_plot_emparejado_", subset_obj$nombre_subset , ".png")
      ggplot2::ggsave(
        filename = file.path(path_results, file_name),
        plot = subset_obj$exp_plot_emparejado,
        width = 15, height = (3 * n_rows) + 2, dpi = 300, scale = 0.8
      )
      message(paste0("Guardado gráficos emparejados con la expresion de genes significativos: ", file_name))
    }
  }
  return(subset_obj)
}




# GUARDAR TABLA DE EXPRESION COMO HOJA DE CALCULO .CSV
guardar_tabla_expresion <- function(subset_obj, codigos_genes_df){
  subset_obj$results_df_sig_exportar <- subset_obj$results_df_sig %>%
    transmute(
      "gen" = gene,
      "media de conteos normalizados" = signif(baseMean, 3),
      "log2 fold change" = signif(log2FoldChange, 3),
      "error estandar de LFC" = signif(lfcSE, 2),
      "estadistico de Wald" = signif(stat, 3),
      # pvalue: redondear a 3 digitos; si es muy pequeño, representar < ...
      "p-value" = ifelse(pvalue <0.001, "< 0.001", 
                         format(round(pvalue, digits = 3), scientific = FALSE)),
      # p-adj: redondear a 4 digitos; si es muy pequeño, representar < ...; añadir estrellas significacion
      "p-value ajustado" = paste0(ifelse(padj <0.0001, 
                                         "< 0.0001", 
                                         ifelse(padj < 0.001, 
                                                format(round(padj, digits = 4), scientific = FALSE),
                                                format(round(padj, digits = 3), scientific = FALSE)
                                                )
                                        ),
                                  " (", sig_stars, ")"
                                  )
      )
  file_name = paste0("resultados_genes_significativos__", subset_obj$nombre_subset , ".csv")
  readr::write_excel_csv(subset_obj$results_df_sig_exportar, file.path(path_results, file_name))
  return(subset_obj)
}




## ------------------------------------------------------------------------
## ////////////////////////////////////////////////////////////////////////
##                ANALISIS DE ENRIQUECIMIENTO FUNCIONAL
## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## ------------------------------------------------------------------------

## Compatible con enrichplot 1.26.6

#/////////////////////////////////////////
#---------- GUARDADO Y ESTILOS -----------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


# GUARDAR UN PLOT (ggplot2::ggsave)
# Crea la carpeta (si no existe) y guarda el objeto ggplot
save_plot <- function(plot, filename, dir_out, width = 9, height = 7, dpi = 300){
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dir_out, filename), plot = plot,
                  width = width, height = height, dpi = dpi)
}

# GUARDAR UNA TABLA (readr::write_csv)
# Crea la carpeta (si no existe) y exporta un data.frame/tibble a CSV
save_tbl <- function(df, filename, dir_out){
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, file.path(dir_out, filename))
}

# AJUSTE PARA DOTPLOTS LARGOS
# p: el dotplot
# wrap_chars: nº de caracteres max por linea
# left_margin_pts: margen izquierdo para que quepan titulos largos
wrap_dot <- function(p, wrap_chars = 45, left_margin_pts = 160) {   
  p +
    # Divide etiquetas del eje Y en varias lineas (str_wrap)
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_chars)) +
    coord_cartesian(clip = "off") +
    scale_color_viridis_c(name = "p.adjust", breaks = scales::pretty_breaks(5)) +   # Color de p.adjust en escala viridis
    theme(plot.margin = margin(10, 20, 10, left_margin_pts),   
          axis.text.y = element_text(size = 11))
}

# CORREGIR NOMBRES DE ARCHIVO
# maxlen: limita longitud maxima del nombre
# Sustituye caracteres problemáticos /\:*?"<>| por _
safe_filename <- function(x, maxlen = 60){
  x <- gsub("[^[:alnum:]_\\-]+", "_", x)
  if (nchar(x) > maxlen) substr(x, 1, maxlen) else x
}

# GUARDAR GSEAPLOT2 DE TOP N TERMINOS
# gsea_obj: procedente de gseGO, gseKEGG, gsePathway
# label_prefix: prefijo de etiqueta de archivo
# top_n: nº de terminos de interes a plotear
save_top_gsea_plots <- function(gsea_obj, label_prefix, out_dir, top_n = 3, width = 8, height = 6){
  # Obtener tabla de resultados y filtrar vacios
  df <- as.data.frame(gsea_obj)
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  # Ordena terminos por FDR y, si hay empate y si existe, tambien por |NES|
  if ("NES" %in% names(df)) {
    df <- df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
  } else {
    df <- df[order(df$p.adjust), , drop = FALSE]
  }
  n <- min(nrow(df), top_n)   # Para no pedir mas graficos de los que hay
  # Generar y guardar cada gseaplot2
  for (i in seq_len(n)) {
    id <- df$ID[i]
    ttl <- df$Description[i]
    p  <- enrichplot::gseaplot2(gsea_obj, geneSetID = id, title = ttl)
    fn <- paste0(label_prefix, "_gseaplot2_TOP", i, "_", safe_filename(ttl), ".png")
    save_plot(p, fn, out_dir, width = width, height = height)
  }
}



#/////////////////////////////////////////
#------- ENRIQUECIMIENTO FUNCIONAL -------
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

# EJECUTAR ANALISIS DE ENRIQUECIMIENTO FUNCIONAL
# Ejecuta GSEA y ORA, frente a bases de datos GO:BP, KEGG, Reactome
# Argumentos:
#   subset_obj: contiene elemento $results_df con columnas: gene, log2FoldChange, pvalue, padj, stat (Wald)
#   label: etiqueta de archivos
#   OrgDb: libreria de anotaciones de Bioconductor -> org.Mm.eg.db (ratón), org.Hs.eg.db (humano)
#   kegg_org: codigo de especie para KEGG -> "mmu" ratón, "hsa" humano
#   reactome_org: especie para Reactome -> "mouse", "humano"; con NULL se infiere de kegg_org
#   out_root: carpeta raiz para guardar archivos de enriquecimiento
#   fdr_thr: umbral de FDR para ORA
#   gsea_min_size, gsea_max_size: rango de genes de los terminos GSEA que nos interesan
#       (terminos muy pequeños son inestables, terminos muy amplios son demasiado genericos)
#   top_n_gsea_terms: nº de terminos para generar curvas de enriquecimiento gseaplot2

run_enrichment <- function(subset_obj,
                           label,
                           OrgDb,                 
                           kegg_org = "mmu",
                           reactome_org = NULL,
                           out_root = file.path(path_results, "Enriquecimiento_Funcional"),
                           fdr_thr = 0.05,
                           gsea_min_size = 10, gsea_max_size = 500,
                           top_n_gsea_terms = 3) {
  
  # 1) PREPARAR DATOS
  
  # Asignar especie para Reactome si NULL
  if (is.null(reactome_org)) {
    reactome_org <- if (tolower(kegg_org) %in% c("hsa","human")) "human" else "mouse"
  }
  
  # Detener analisis si no existe $results_df
  stopifnot(!is.null(subset_obj$results_df)) 
  # Cargar $results_df en objeto res_all y asignar columna de gene si genes en nombres de fila
  res_all <- subset_obj$results_df
  if (is.null(res_all$gene)) res_all$gene <- rownames(res_all)
  
  # Asegura columna 'stat' para ranking GSEA; si falta, crea proxy usando pvalue y signo de LFC
  if (!("stat" %in% names(res_all)) || all(is.na(res_all$stat))) {
    if (!("pvalue" %in% names(res_all))) stop("No hay 'stat' ni 'pvalue' en results_df para construir ranking.")
    res_all$stat <- with(res_all, sign(log2FoldChange) * -log10(pmax(pvalue, .Machine$double.xmin)))  # evita -Inf
  }
  
  # Directorio de salida
  out_dir <- file.path(out_root, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Ranking para GSEA (estadistico de Wald o proxy); descartar si NA
  rank_df <- res_all |>
    dplyr::select(gene, stat) |>
    dplyr::filter(!is.na(gene), !is.na(stat))
  
  # Conversion de genes SYMBOL -> ENTREZID
  conv_all <- bitr(rank_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Resolver duplicados de ENTREZID (conservar el de mayor |stat|)
  rank_merged <- rank_df |>
    dplyr::inner_join(conv_all, by = c("gene" = "SYMBOL")) |>
    dplyr::group_by(ENTREZID) |>
    dplyr::slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  # Contruir geneList para gse*: vector nombrado por ENTREZID, en orden descendente
  geneList <- rank_merged$stat
  names(geneList) <- rank_merged$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Definir universo de genes para ORA
  universe_entrez <- unique(conv_all$ENTREZID)
  
  
  # -----  
  # 2) GSEA (GO:BP, KEGG, Reactome)

  # Aplicar GSEA sobre GO:BP
  gse_go_bp <- gseGO(
    geneList      = geneList,      # Lista de genes
    OrgDb         = OrgDb,         # Libreria de anotaciones (especie)
    keyType       = "ENTREZID",    # Clave acorde a geneList
    ont           = "BP",          # Ontologia dentro de GO
    minGSSize     = gsea_min_size, # Filtros para tamaño del conjunto
    maxGSSize     = gsea_max_size, # "
    pAdjustMethod = "BH",          # Metodo de correccion por tests multiples (FDR)
    verbose       = FALSE          # Mensajes para consola
  )

  # Aplicar GSEA sobre KEGG
  gse_kegg <- gseKEGG(
    geneList      = geneList,      # Lista de genes
    organism      = kegg_org,      # Especie
    minGSSize     = gsea_min_size, # Filtro para tamaño del conjunto
    pAdjustMethod = "BH",          # Metodo de correccion por tests multiples (FDR)
    verbose       = FALSE          # Mensajes para consola
  )
  
  # Aplicar GSEA sobre Reactome
  gse_react <- ReactomePA::gsePathway(
    geneList      = geneList,      # Lista de genes
    organism      = reactome_org,  # Especie
    minGSSize     = gsea_min_size, # Filtro para tamaño del conjunto
    pAdjustMethod = "BH",          # Metodo de correccion por tests multiples (FDR)
    verbose       = FALSE          # Mensajes para consola
  )
  
  # Exportar tablas GSEA con resultados
  if (!is.null(gse_go_bp)  && nrow(as.data.frame(gse_go_bp))  > 0) save_tbl(as.data.frame(gse_go_bp),  paste0("GSEA_GO_BP_",    label, ".csv"), out_dir)
  if (!is.null(gse_kegg)   && nrow(as.data.frame(gse_kegg))   > 0) save_tbl(as.data.frame(gse_kegg),   paste0("GSEA_KEGG_",     label, ".csv"), out_dir)
  if (!is.null(gse_react)  && nrow(as.data.frame(gse_react))  > 0) save_tbl(as.data.frame(gse_react),  paste0("GSEA_Reactome_", label, ".csv"), out_dir)
  
  
  # -----  
  # 3) ORA (GO:BP, KEGG, Reactome) UP/DOWN
  
  # Separar genes UP y DOWN a partir de resultados ajustados
  genes_up   <- res_all |> dplyr::filter(!is.na(padj), padj < fdr_thr, log2FoldChange > 0) |> dplyr::pull(gene)
  genes_down <- res_all |> dplyr::filter(!is.na(padj), padj < fdr_thr, log2FoldChange < 0) |> dplyr::pull(gene)

  # Aplicar ORA sobre GO:BP de genes UP y DOWN
  ego_up <- enrichGO(gene = genes_up, universe = rank_df$gene,
                     OrgDb = OrgDb, keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = "BH", qvalueCutoff = fdr_thr)
  ego_down <- enrichGO(gene = genes_down, universe = rank_df$gene,
                       OrgDb = OrgDb, keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = "BH", qvalueCutoff = fdr_thr)
  
  # Simplificar redundancia semantica en GO
  if (!is.null(ego_up)   && nrow(as.data.frame(ego_up))   > 0)   ego_up   <- simplify(ego_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
  if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0)   ego_down <- simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # Mapear genes UP y DOWN a ENTREZ, para KEGG y Reactome
  map_up   <- if (length(genes_up))   bitr(genes_up,   fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb) else NULL
  map_down <- if (length(genes_down)) bitr(genes_down, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb) else NULL
  
  # Aplicar ORA sobre KEGG de genes UP y DOWN
  ekegg_up <- if (!is.null(map_up) && nrow(map_up) > 0) {
    enrichKEGG(gene = unique(map_up$ENTREZID), universe = universe_entrez, organism = kegg_org, pvalueCutoff = 0.05)
  } else NULL
  ekegg_down <- if (!is.null(map_down) && nrow(map_down) > 0) {
    enrichKEGG(gene = unique(map_down$ENTREZID), universe = universe_entrez, organism = kegg_org, pvalueCutoff = 0.05)
  } else NULL
  
  # Aplicar ORA sobre Reactome de genes UP y DOWN
  ereact_up <- if (!is.null(map_up) && nrow(map_up) > 0) {
    ReactomePA::enrichPathway(gene = unique(map_up$ENTREZID), universe = universe_entrez,
                              organism = reactome_org, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  } else NULL
  ereact_down <- if (!is.null(map_down) && nrow(map_down) > 0) {
    ReactomePA::enrichPathway(gene = unique(map_down$ENTREZID), universe = universe_entrez,
                              organism = reactome_org, pvalueCutoff = 0.05, pAdjustMethod = "BH")
  } else NULL
  
  # Exportar tablas ORA (si existen)
  if (!is.null(ego_up)      && nrow(as.data.frame(ego_up))      > 0) save_tbl(as.data.frame(ego_up),      paste0("ORA_GO_BP_UP_",      label, ".csv"), out_dir)
  if (!is.null(ego_down)    && nrow(as.data.frame(ego_down))    > 0) save_tbl(as.data.frame(ego_down),    paste0("ORA_GO_BP_DOWN_",    label, ".csv"), out_dir)
  if (!is.null(ekegg_up)    && nrow(as.data.frame(ekegg_up))    > 0) save_tbl(as.data.frame(ekegg_up),    paste0("ORA_KEGG_UP_",       label, ".csv"), out_dir)
  if (!is.null(ekegg_down)  && nrow(as.data.frame(ekegg_down))  > 0) save_tbl(as.data.frame(ekegg_down),  paste0("ORA_KEGG_DOWN_",     label, ".csv"), out_dir)
  if (!is.null(ereact_up)   && nrow(as.data.frame(ereact_up))   > 0) save_tbl(as.data.frame(ereact_up),   paste0("ORA_Reactome_UP_",   label, ".csv"), out_dir)
  if (!is.null(ereact_down) && nrow(as.data.frame(ereact_down)) > 0) save_tbl(as.data.frame(ereact_down), paste0("ORA_Reactome_DOWN_", label, ".csv"), out_dir)
  

  # -----  
  # 3) PLOTS
  #    Para GSEA: dotplot, ridgeplot, GSEAplot2 (top N terminos)
  #    Para ORA: dotplot, emapplot (solo GO) y cnetplot (solo GO)
  
  # Plots GSEA GO
  if (!is.null(gse_go_bp) && nrow(as.data.frame(gse_go_bp)) > 0) {
    # Dotplot    
    p1 <- enrichplot::dotplot(gse_go_bp, x = "NES", showCategory = 20) +
      ggtitle(paste0("GSEA — GO BP — ", label))
    save_plot(wrap_dot(p1), paste0("GSEA_GO_BP_dotplot_", label, ".png"), out_dir)
    # Ridgeplot
    p2 <- enrichplot::ridgeplot(gse_go_bp, showCategory = 15) +
      scale_fill_viridis_c(name = "p.adjust", breaks = scales::pretty_breaks(5)) +
      ggtitle(paste0("GSEA — GO BP (ridge) — ", label))
    save_plot(p2, paste0("GSEA_GO_BP_ridgeplot_", label, ".png"), out_dir, width = 10, height = 7)
    # GSEAplot2 (top N terminos)
    save_top_gsea_plots(gse_go_bp, paste0("GSEA_GO_BP_", label), out_dir, top_n = top_n_gsea_terms)
  }
  
  # Plots GSEA KEGG
  if (!is.null(gse_kegg) && nrow(as.data.frame(gse_kegg)) > 0) {
    # Dotplot
    p4 <- enrichplot::dotplot(gse_kegg, x = "NES", showCategory = 20) +
      ggtitle(paste0("GSEA — KEGG — ", label))
    save_plot(wrap_dot(p4), paste0("GSEA_KEGG_dotplot_", label, ".png"), out_dir)
    # Ridgeplot
    p_rk <- enrichplot::ridgeplot(gse_kegg, showCategory = 15) +
      scale_fill_viridis_c(name = "p.adjust", breaks = scales::pretty_breaks(5)) +
      ggtitle(paste0("GSEA — KEGG (ridge) — ", label))
    save_plot(p_rk, paste0("GSEA_KEGG_ridgeplot_", label, ".png"), out_dir, width = 10, height = 7)
    # GSEAplot2 (top N terminos)
    save_top_gsea_plots(gse_kegg, paste0("GSEA_KEGG_", label), out_dir, top_n = top_n_gsea_terms)
  }
  
  # Plots GSEA Reactome
  if (!is.null(gse_react) && nrow(as.data.frame(gse_react)) > 0) {
    # Dotplot
    p6 <- enrichplot::dotplot(gse_react, x = "NES", showCategory = 20) +
      ggtitle(paste0("GSEA — Reactome — ", label))
    save_plot(wrap_dot(p6), paste0("GSEA_Reactome_dotplot_", label, ".png"), out_dir)
    # Ridgeplot
    p_rr <- enrichplot::ridgeplot(gse_react, showCategory = 15) +
      scale_fill_viridis_c(name = "p.adjust", breaks = scales::pretty_breaks(5)) +
      ggtitle(paste0("GSEA — Reactome (ridge) — ", label))
    save_plot(p_rr, paste0("GSEA_Reactome_ridgeplot_", label, ".png"), out_dir, width = 10, height = 7)
    # GSEAplot2 (top N terminos)
    save_top_gsea_plots(gse_react, paste0("GSEA_Reactome_", label), out_dir, top_n = top_n_gsea_terms)
  }
  
  # Plots ORA GO (UP) + redes
  if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
    # Dotplot
    p_up <- enrichplot::dotplot(ego_up, showCategory = 15) + ggtitle(paste0("ORA — GO BP (UP) — ", label))
    save_plot(wrap_dot(p_up), paste0("ORA_GO_BP_UP_dotplot_", label, ".png"), out_dir)
    # emapplot
    ego_up_sim <- pairwise_termsim(ego_up)
    p_emap <- enrichplot::emapplot(ego_up_sim, showCategory = 30, layout = "nicely", label_format = 40) +
      ggraph::scale_edge_alpha(range = c(0.2, 0.8)) +
      scale_color_viridis_c(end = 0.9) +
      ggtitle(paste0("Enrichment map — GO BP UP — ", label))
    save_plot(p_emap, paste0("ORA_GO_BP_UP_emapplot_", label, ".png"), out_dir, width = 11, height = 8)
    # cnetplot
    fc <- setNames(res_all$log2FoldChange, res_all$gene)
    p_cnet <- enrichplot::cnetplot(ego_up, showCategory = 10, foldChange = fc,
                                   circular = FALSE, layout = "nicely",
                                   node_label = "all", label_format = 35) +
      scale_color_viridis_c(end = 0.9) +
      ggtitle(paste0("Category netplot — GO BP UP — ", label))
    save_plot(p_cnet, paste0("ORA_GO_BP_UP_cnetplot_", label, ".png"), out_dir, width = 11, height = 8)
  }
  
  # Plots ORA GO (DOWN) + redes
  if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
    # Dotplot
    p_dn <- enrichplot::dotplot(ego_down, showCategory = 15) + ggtitle(paste0("ORA — GO BP (DOWN) — ", label))
    save_plot(wrap_dot(p_dn), paste0("ORA_GO_BP_DOWN_dotplot_", label, ".png"), out_dir)
    # emapplot
    ego_dn_sim <- pairwise_termsim(ego_down)
    p_emap_dn <- enrichplot::emapplot(ego_dn_sim, showCategory = 30, layout = "nicely", label_format = 40) +
      ggraph::scale_edge_alpha(range = c(0.2, 0.8)) +
      scale_color_viridis_c(end = 0.9) +
      ggtitle(paste0("Enrichment map — GO BP DOWN — ", label))
    save_plot(p_emap_dn, paste0("ORA_GO_BP_DOWN_emapplot_", label, ".png"), out_dir, width = 11, height = 8)
    # cnetplot
    fc <- setNames(res_all$log2FoldChange, res_all$gene)
    p_cnet_dn <- enrichplot::cnetplot(ego_down, showCategory = 10, foldChange = fc,
                                      circular = FALSE, layout = "nicely",
                                      node_label = "all", label_format = 35) +
      scale_color_viridis_c(end = 0.9) +
      ggtitle(paste0("Category netplot — GO BP DOWN — ", label))
    save_plot(p_cnet_dn, paste0("ORA_GO_BP_DOWN_cnetplot_", label, ".png"), out_dir, width = 11, height = 8)
  }
  
  # Plots ORA KEGG
  # Dotplot UP
  if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) {
    p_kup <- enrichplot::dotplot(ekegg_up, showCategory = 15) + ggtitle(paste0("ORA — KEGG (UP) — ", label))
    save_plot(wrap_dot(p_kup), paste0("ORA_KEGG_UP_dotplot_", label, ".png"), out_dir)
  }
  # Dotplot DOWN
  if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) {
    p_kdn <- enrichplot::dotplot(ekegg_down, showCategory = 15) + ggtitle(paste0("ORA — KEGG (DOWN) — ", label))
    save_plot(wrap_dot(p_kdn), paste0("ORA_KEGG_DOWN_dotplot_", label, ".png"), out_dir)
  }
  
  # Plots ORA Reactome
  # Dotplot UP
  if (!is.null(ereact_up) && nrow(as.data.frame(ereact_up)) > 0) {
    p_ru <- enrichplot::dotplot(ereact_up, showCategory = 15) + ggtitle(paste0("ORA — Reactome (UP) — ", label))
    save_plot(wrap_dot(p_ru), paste0("ORA_Reactome_UP_dotplot_", label, ".png"), out_dir)
  }
  # Dotplot DOWN
  if (!is.null(ereact_down) && nrow(as.data.frame(ereact_down)) > 0) {
    p_rd <- enrichplot::dotplot(ereact_down, showCategory = 15) + ggtitle(paste0("ORA — Reactome (DOWN) — ", label))
    save_plot(wrap_dot(p_rd), paste0("ORA_Reactome_DOWN_dotplot_", label, ".png"), out_dir)
  }
  
  # Salida
  message(">> Enriquecimiento completado: ", out_dir)
  invisible(list(
    gse_go_bp  = gse_go_bp,
    gse_kegg   = gse_kegg,
    gse_react  = gse_react,
    ego_up     = ego_up,
    ego_down   = ego_down,
    ekegg_up   = ekegg_up,
    ekegg_down = ekegg_down,
    ereact_up  = ereact_up,
    ereact_down= ereact_down
  ))
}

