# TFM_bioinformatica
Trabajo Fin de Máster, Máster en Bioinformática  
Universidad Internacional de la Rioja (UNIR)  
Autor: Víctor Saavedra Yturriagagoitia  
Fecha: Septiembre 2025  


Este repositorio acompaña el Trabajo de Fin de Máster *“Análisis transcriptómico y caracterización funcional de biomarcadores inflamatorios en respuesta a la suplementación con ácidos grasos omega-3”*. Reúne el código en R (funciones y pipelines), los datasets de RNA-seq descargados de GEO ([GSE132648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132648) para humano; [GSE153648](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153648) para ratón) y todas las figuras y tablas generadas.  
  
El flujo en R reproduce: preparación de datos y metadatos, DESeq2 para expresión diferencial, visualización (correlación, PCA, volcano, heatmaps, plots de expresión) y enriquecimiento funcional (GO/KEGG/Reactome) con clusterProfiler, ReactomePA y enrichplot.  

> For a quick visual overview of the main RNA-seq analyses and figures, see [RNAseq_figure_showcase.md](RNAseq_figure_showcase.md).

## Índice
- [Estructura del repositorio](#estructura-del-repositorio)
- [Requisitos](#requisitos)
- [Flujo del análisis](#flujo-del-análisis)
- [Uso rápido](#uso-rápido)
- [Datasets](#datasets)
- [Results](#results)
- [Scripts](#scripts)
- [Citación](#citación)



## Estructura del repositorio

| Carpeta  | Contenido                       |
|--------- |---------------------------------|
| **Datasets/** | Datos de RNA-seq descargados de GEO (conteos sin normalizar). |
| **Results/**  | Figuras y tablas generadas por los scripts. |
| **Scripts/**  | Código en R creado en este trabajo (funciones y pipelines para cada dataset). |
| **FIGURAS_EJEMPLO**  | Algunas figuras de ejemplo para visualización rápida. |
| **README.md** | Este documento. |
| **TFM_VictorSaavedra.pdf** | Documento del Trabajo de Fin de Máster. |


### Vista en árbol

La nomenclatura de las figuras se detalla en la sección [Results](#results).

```
TFM_bioinformatica/
├─ Datasets/
│  ├─ GSE132648_Souza2019/
│  │  ├─ <conteos>.txt
│  │  └─ <diseño>.png
│  └─ GSE153648_Sorokin2023/
│     ├─ <conteos>.xlsx
│     └─ <diseño>.png
├─ Results/
│  ├─ Results_GSE132648_Souza2019/
│  │  ├─ resultados_genes_significativos__all.csv
│  │  ├─ (*) <figuras_deseq2>.png
│  │  └─ Enriquecimiento_Funcional/
│  │     └─ Souza_aceitevsplacebo/
│  │        ├─ (*) <tablas_GSEA/ORA>.csv
│  │        └─ (*) <figuras_GSEA/ORA>.png
│  ├─ Results_GSE153648_Sorokin2023/
│  │  ├─ (*) resultados_genes_significativos__<tejido>_<trat1vstrat2>.csv
│  │  ├─ (*) <figuras_deseq2>.png
│  │  ├─ venn_aorta_genessig.png
│  │  └─ Enriquecimiento_Funcional/
│  │     └─ (*) subset_<tejido>_<trat1vstrat2>/
│  │        ├─ (*) <tablas_GSEA/ORA>.csv
│  │        └─ (*) <figuras_GSEA/ORA>.png
│  └─ resultados_genes_significativos.xlsx
├─ Scripts/
│  ├─ RNAseq_GSE132648_Souza2019.R
│  ├─ RNAseq_GSE153648_Sorokin2023.R
│  └─ RNAseq_funciones.R
├─ FIGURAS_EJEMPLO.md
├─ README.md
└─ TFM_VictorSaavedra.pdf
```

Nota: las ramas con (*) indican multiplicidad de archivos o carpetas:

<br>  

`(*) <figuras_deseq2>.png`  
-Correlación entre muestras:  mapa de calor de correlación,  PCA plot,  PCA 3D plot.  
-Dispersión del conjunto de datos:  plot de estimaciones de dispersión,  plot de media vs varianza.  
-Resultados de expresión diferencial:  mapa de calor de expresión,  plot de expresión de top N genes,  plot de expresión emparejado de top N genes,  volcano plot,  volcano plot con etiquetas de genes significativos.  
<details>
    <summary>(ver detalle de figuras) </summary>
    
```
├─ correlation_heatmap_*.png
├─ estimateddispersion_plot_*.png
├─ expression_heatmap_*.png
├─ expression_plot_*_<top_n_genes>.png
├─ expression_plot_emparejado_*.png
├─ meanvar_plot_*.png
├─ pca3D_plot_*.png
├─ pca_plot_*.png
├─ volcano_plot_*.png
└─ volcano_plot_labels_*.png

* distinto para cada subset (tejido y combinación de tratamientos)  
```
</details><br>

`(*) <tablas_GSEA/ORA>.csv`  
-GSEA para GO:BP, KEGG, Reactome;  
-ORA de genes DOWN y UP para GO:BP, KEGG, Reactome
<details>
    <summary>(ver detalle de tablas) </summary>
    
```
├─ GSEA_GO_BP_*.csv
├─ GSEA_KEGG_*.csv
├─ GSEA_Reactome_*.csv
├─ ORA_GO_BP_DOWN_*.csv
├─ ORA_GO_BP_UP_*.csv
├─ ORA_KEGG_DOWN_*.csv
├─ ORA_KEGG_UP_*.csv
├─ ORA_Reactome_DOWN_*.csv
└─ ORA_Reactome_UP_*.csv

* distinto para cada subset (tejido y combinación de tratamientos)  
```
</details><br>

`(*) <figuras_GSEA/ORA>.png`  
-GSEA para GO:BP, KEGG, Reactome: dotplot, ridgeplot, gseaplot2 de N top términos;   
-ORA de genes DOWN y UP para GO:BP, KEGG, Reactome: dotplot, cnetplot(solo GO:BP), emapplot (solo GO:BP)
<details>
    <summary>(ver detalle de figuras) </summary>
    
```
├─ GSEA_GO_BP_dotplot_*.png
├─ GSEA_GO_BP_ridgeplot_*.png
├─ (*) GSEA_GO_BP_*_gseaplot2_TOP#_<termino>.png
├─ GSEA_KEGG_dotplot_*.png
├─ GSEA_KEGG_ridgeplot_*.png
├─ (*) GSEA_KEGG_*_gseaplot2_TOP#_<termino>.png
├─ GSEA_Reactome_dotplot_*.png
├─ GSEA_Reactome_ridgeplot_*.png
├─ (*) GSEA_Reactome_*_gseaplot2_TOP#_<termino>.png
├─ ORA_GO_BP_DOWN_cnetplot_*.png
├─ ORA_GO_BP_DOWN_dotplot_*.png
├─ ORA_GO_BP_DOWN_emapplot_*.png
├─ ORA_GO_BP_UP_cnetplot_*.png
├─ ORA_GO_BP_UP_dotplot_*.png
├─ ORA_GO_BP_UP_emapplot_*.png
├─ ORA_KEGG_DOWN_dotplot_*.png
├─ ORA_KEGG_UP_dotplot_*.png
├─ ORA_Reactome_DOWN_dotplot_*.png
└─ ORA_Reactome_UP_dotplot_*.png

* distinto para cada subset (tejido y combinación de tratamientos)
```
</details><br>

Para más información sobre los *subsets* (tejido y combinación de tratamiento), ver la **Arquitectura del código (*subsets* como listas)** en la sección [Scripts](#scripts).

---

## Requisitos

R (versión 4.4.2)  
BiocManager (3.20)  
Bioconductor (3.20)  

**Paquetes Bioconductor:**  
DESeq2 (1.46.0), GEOquery (2.74.0), SummarizedExperiment (1.36.0), biomaRt (2.62.1), org.Mm.eg.db (3.20.0), org.Hs.eg.db (3.20.0), clusterProfiler (4.14.6), enrichplot (1.26.6), ReactomePA (1.50.0).  

**Paquetes CRAN:**  
pheatmap (1.0.12), ggplot2 (3.5.2), dplyr (1.1.4), tidyr (1.3.1), tibble (3.2.1), readr (2.1.5), stringr (1.5.1), ggrepel (0.9.6), patchwork (1.3.0), plotly (4.10.4), scales (1.3.0), ggraph (2.2.2) y ggridges (0.5.6).  



---

## Flujo del análisis


```
[Conteos crudos de GEO]
[Metadatos de GEO]
   ↓ GESTION SUBSETS Y METADATOS
[Datos y metadatos adaptados]
   | ANÁLISIS DE EXPRESIÓN DIFERENCIAL con DESeq2
   | (EXPLORACIÓN DE LOS DATOS)
   ├─→ [Figuras: Correlación entre muestras] ─→ heatmap correlacion, PCA, PCA 3D
   | (MODELAJE DE EXPRESIÓN DIFERENCIAL)
   ├─→ [Figuras: Dispersión de datos] ─→ estimaciones de dispersión, media vs varianza
   ↓ (CÁLCULO DE RESULTADOS)
A. [Métricas estadísticas para cada gen (baseMean, log2FC, lfcSE, stat Wald, pvalue, padj)]
B. [Lista de genes significativos UP y DOWN]
C. [Conteos normalizados]
   |   | VISUALIZACIÓN DE RESULTADOS
   |   ├─→ [Tabla: Genes significativos con métricas estadísticas]
   |   └─→ [Figuras: Resultados de expr. diferencial] ─→ heatmap, volcano plot, plot de expresión, plot de expresión emparejado
   | ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL
   A ─→ |ranking stat| ─→ [Figuras y tablas GSEA: GO/KEGG/Reactome]
   B ─→ |UP/DOWN| ─→ [Figuras y tablas ORA: GO/KEGG/Reactome]
```

---
## Uso rápido

1. Clonar el repositorio:  
   ```bash
   git clone https://github.com/vsaavedra927/TFM_bioinformatica.git
   cd TFM_bioinformatica
   
2. Instalar paquetes (ver sección [Requisitos](#requisitos)).  
Nota: Este paso es opcional: los paquetes se instalan automáticamente al ejecutar el pipeline de funciones en el paso 5.  
Desde R:
    ```r
    source("Scripts/RNAseq_funciones.R"); load_pkgs() 
  
4. Colocar los datasets en `Datasets/` (ya incluidos).  
  
5. Ejecutar los pipelines de la carpeta `Scripts/`. Los resultados se guardan automáticamente en `Results/*` (figuras y CSV).  
Desde R:
    ```r
    source("Scripts/RNAseq_funciones.R")
    source("Scripts/RNAseq_GSE132648_Souza2019.R")
    source("Scripts/RNAseq_GSE153648_Sorokin2023.R")
  
6. Variables y parámetros útiles para configurar:  
Ruta de carpeta Github: `path_carpetaGithub`  
Umbrales FDR y |log2FC|: `umbral_fdr`, `umbral_lfc`  
Nota: Las funciones incluyen numerosos parámetros para analizar con distintos tests estadísticos y para generar figuras de distintos formatos. Para más información, consultar los comentarios de los scripts y la documentación de las librerías.  


---

    
## Datasets
Entrada (input) de los scripts.  
Contiene los conteos crudos de RNA-seq sin editar, extraídos de GEO, y un esquema del diseño experimental en ambos estudios.  
  
Corresponden a los siguientes conjuntos de datos de GEO:  
- **GSE132648, publicado por Souza et al. (2019):** voluntarios humanos sanos suplementados con aceite marino enriquecido en omega-3 o placebo; muestras de sangre periférica.  
- **GSE153648, publicado por Sorokin et al. (2023):** ratones suplementados con DHA, EPA o placebo; muestras de aorta, hígado, piel.  



## Results 
Salida (output) de los scripts.  
Contiene las figuras y tablas resultantes de los análisis de expresión diferencial y de enriquecimiento funcional.
Contenido de la carpeta `Results/`:
- `Results_GSE132648_Souza2019/`: carpeta con todos los resultados sobre los datos de humano.  
- `Results_GSE153648_Sorokin2023/`: carpeta con todos los resultados sobre los datos de ratón.  
- `resultados_genes_significativos_SorokinSouza.xlsx`: archivo Excel recopilatorio con todos los genes significativos del estudio.    
     
Nomenclatura de archivos de resultados:
-  Figuras de correlación entre muestras (`.png`) (demuestran homogeneidad de las muestras):  
    -  heatmap de correlacion: `correlation_heatmap_*.png`   
    -  PCA: `pca_plot_*.png`  
    -  PCA 3D: `pca3D_plot_*.png`  
-  Figuras de dispersión de datos (`.png`) (demuestran aptitud del modelo estadístico):  
    -  estimaciones de dispersión: `estimateddispersion_plot_*.png`  
    -  media vs varianza: `meanvar_plot_*.png`  
-  Figuras de resultados de expresión diferencial (`.png`):  
    -  mapa de calor de expresión: `expression_heatmap_*.png`. En alta resolución con sufijo `_highres.png`.
    -  plot de expresión de top N genes: `expression_plot_*_<topNgenes>.png`   
    -  plot de expresión emparejado de top N genes: `expression_plot_emparejado_*.png`   
    -  volcano plot: `volcano_plot_*.png`  
    -  volcano plot con etiquetas de genes significativos: `volcano_plot_labels_*.png`
-  Tabla de genes significativos con métricas estadísticas: `resultados_genes_significativos__*.csv`  
-  Carpeta `Enriquecimiento_Funcional/`: resultados de enriquecimiento (figuras y tablas; GSEA y ORA para GO:BP, KEGG, Reactome)   
    -  Figuras(`.png`) de Gene Set Enrichment Analysis (GSEA) para GO:BP, KEGG, Reactome; para cada base de datos se incluyen dotplot, ridgeplot y gseaplot2 de top N términos:  
        - `GSEA_<GO_BP/KEGG/Reactome>_dotplot_*.png`  
        - `GSEA_<GO_BP/KEGG/Reactome>_ridgeplot_*.png`  
        - `GSEA_<GO_BP/KEGG/Reactome>_*_gseaplot2_TOP<#>_<termino>.png`  
    -  Figuras (`.png`) de Over-Representation Analysis (ORA) de genes UP y DOWN para GO:BP, KEGG, Reactome; para cada base de datos se incluyen dotplot; para GO:BP además cnetplot y emapplot:  
        - `ORA_<GO_BP/KEGG/Reactome>_<UP/DOWN>_dotplot_*.png`  
        - `ORA_GO_BP_<UP/DOWN>_cnetplot_*.png`  
        - `ORA_GO_BP_<UP/DOWN>_emapplot_*.png`  

**Nota:** los `*` en la lista de nombres de archivo indica multiplicidad de tejidos y comparaciones:  
Los resultados de Souza son para un solo tejido (sangre periférica de humano) y una sola comparación (aceite vs placebo).  
Los resultados de Sorokin son para tres tejidos (aorta, hígado y piel de ratón) y tres comparaciones (DHA vs ctrl, EPA vs ctrl, EPA vs DHA).  



## Scripts
Scripts de R.  
Estos scripts usan como entrada (input) los archivos de la carpeta Datasets y generan como salida (output) los archivos de la carpeta Results.  
Contenido de la carpeta `Scripts/`:
- `RNAseq_funciones.R`: script con las funciones a utilizar.  
- `RNAseq_GSE132648_Souza2019.R`: script para analizar datos de Souza et al. (humano).  
- `RNAseq_GSE153648_Sorokin2023.R`: script para analizar datos de Sorokin et al. (ratón).  


  
### Arquitectura del código (*subsets* como listas)

Cada *subset* (uno por tejido y comparación) se implementa como una lista de R, que incluye los datos de entrada (conteos y metadatos), información relevante sobre ese *subset*, y todos los objetos generados durante el análisis. Las funciones reciben un `subset_obj` (la lista correspondiente a un *subset*), trabajan sobre él y devuelven el mismo `subset_obj` actualizado (incluyendo nuevos elementos). De este modo, el estado queda encapsulado y se asegura la reproducibilildad del flujo de trabajo para los distintos *subsets*.  

Campos típicos de un `subset_obj` (lista no exhaustiva):  
-Información base: `nombre_subset`, `gse_nr`, `especie`, `tejido`, `tratamiento`, `muestras`, `metadata_df`, `count_matrix`.  
-Preprocesado/Exploración: `dds`, `count_matrix_norm`, `vsd`, `vsd_mat`, `vsd_cor`, `meanvar_df`, `pca_plot`.  
-Análisis de expresión diferencial: `results` / `results_lfcshrink`, `results_df`, `results_df_sig`, `count_matrix_norm_sig`.  
-Figuras: `volcano_plot`, `exp_plot`, `exp_plot_emparejado`.  

Patrón de uso (ejemplo):  
  ```r
  subset_obj <- crear_DESeqDataSet(subset_obj, variable_diseno = "tratamiento")
  subset_obj <- normalizar_conteos(subset_obj)
  subset_obj <- estudiar_correlacion_datos(subset_obj, heatmap = TRUE, pca = TRUE, colores = set_colores)
  
  subset_obj <- aplicar_deseq(subset_obj)
  subset_obj <- extraer_resultados_log2fcshrink(
    subset_obj,
    coeficiente = "tratamiento_DHAsup_vs_Control",
    metadatos_global = metadata_df_all,
    reducir_lfc = FALSE, umbral_fdr = 0.05, umbral_lfc = 0.585
  )
  
  subset_obj <- visualizar_resultados(
    subset_obj, metadatos_global = metadata_df_all, col_metadata = col_metadata_para_exp_plot,
    umbral_fdr = 0.05, umbral_lfc = 0.585, colores = set_colores,
    heatmap = TRUE, volcano_plot = TRUE, exp_plot = TRUE, exp_plot_emparejado = FALSE
  )
  ```

Ventajas de este diseño:  
-Reproducibilidad: cada paso deja sus resultados dentro del mismo objeto.  
-Modularidad: se pueden activar/desactivar bloques sin romper el resto.  
-Inspección sencilla: `str(subset_obj)` muestra todo el estado del análisis de ese *subset*; `subject_obj$elemento` muestra un elemento del *subset*.
-Evita efectos colaterales: se reasigna explícitamente la lista devuelta por cada función.  
  
**Nota:** para hacer el código más eficiente, se pueden almacenar todos los nombres de las listas en un vector y acceder a las listas a través de sus nombres en formato *string*; así en los *loops* del pipeline se sigue siempre el patrón leer → actualizar → reasignar el `subset_obj` a partir de sus nombres en formato *string*, para todos los *subsets* por igual.  


---

## Citación

Si usas este repositorio, por favor cita:  
  
**Este repositorio**  
  Saavedra Yturriagagoitia, V. (2025). *TFM_bioinformatica* (v1.0). GitHub.  
  Disponible en: https://github.com/vsaavedra927/TFM_bioinformatica
  
**Conjuntos de datos de GEO**  
- NCBI Gene Expression Omnibus (GEO). GSE132648 — Enriched marine oil supplements increase peripheral blood SPM concentrations and reprogram host immune responses: A randomized double-blind placebo-controlled study. 2019 (actualizado 2020-03-13). Disponible en: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132648. Accedido el: 10-06-2025.  
- NCBI Gene Expression Omnibus (GEO). GSE153648 — Comparison of the dietary omega-3 fatty acids impact on murine psoriasis-like skin inflammation and associated lipid dysfunction. Público el 13-abr-2023. Disponible en: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153648. Accedido el: 20-07-2025.  
  
**Artículos de los autores**  
  - Souza PR, Marques RM, Gomez EA, Colas RA, De Matteis R, Zak A, et al. Enriched Marine Oil Supplements Increase Peripheral Blood Specialized Pro-Resolving Mediators Concentrations and Reprogram Host Immune Responses: A Randomized Double-Blind Placebo-Controlled Study. Circ Res [Internet]. 2020 Jan 3 [cited 2025 May 16];126(1):75–90. Available from: [/doi/pdf/10.1161/CIRCRESAHA.119.315506](https://pubmed.ncbi.nlm.nih.gov/31829100/)
- Sorokin A V., Arnardottir H, Svirydava M, Ng Q, Baumer Y, Berg A, et al. Comparison of  the dietary omega-3 fatty acids impact on murine psoriasis-like skin inflammation and associated lipid dysfunction. Journal of Nutritional Biochemistry [Internet]. 2023 Jul 1 [cited 2025 Jul 23];117. Available from: [https://pubmed.ncbi.nlm.nih.gov/37044136/  ](https://pubmed.ncbi.nlm.nih.gov/37044136/)

