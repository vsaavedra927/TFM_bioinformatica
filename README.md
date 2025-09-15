# TFM_bioinformatica
Trabajo Fin de Máster del Máster en Bioinformática, Universidad Internacional de la Rioja (UNIR)
Autor: Víctor Saavedra Yturriagagoitia
Fecha: Septiembre 2025

Este repositorio contiene conjuntos de RNA-seq obtenidos de Gene Expression Omnibus (GEO), scripts en R para análisis de expresión diferencial (DESeq2) y enriquecimiento funcional (GO/KEGG/Reactome), y los resultados generados.



## Índice
- [Estructura del repositorio](#estructura-del-repositorio)
- [Requisitos](#requisitos)
- [Flujo del análisis](#flujo-del-analisis)
- [Uso rápido](#uso-rápido)
- [Datasets](#datasets)
- [Resultados](#results)
- [Scripts](#scripts)
- [Citación y licencia](#citación-y-licencia)



## Estructura del repositorio

| Carpeta  | Contenido                       |
|--------- |---------------------------------|
| **Datasets/** | Datos de RNA-seq descargados de GEO (conteos sin normalizar). |
| **Results/**  | Figuras y tablas generadas por los scripts. |
| **Scripts/**  | Código en R creado en este trabajo (funciones y pipelines para cada dataset). |
| **README.md** | Este documento. |


### Vista en árbol

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
│  │  └─ Enriquecimiento_Funcional/
│  │     └─ (*) subset_<tejido>_<trat1vstrat2>/
│  │        ├─ (*) <tablas_GSEA/ORA>.csv
│  │        └─ (*) <figuras_GSEA/ORA>.png
│  └─ resultados_genes_significativos.xlsx
├─ Scripts/
│  ├─ RNAseq_GSE132648_Souza2019.R
│  ├─ RNAseq_GSE153648_Sorokin2023.R
│  └─ RNAseq_funciones.R
└─ README.md

Nota: las ramas con (*) indican multiplicidad de archivos o carpetas.  
```


**(*) <figuras_deseq2>.png**  
Correlación entre muestras:  mapa de calor de correlación,  PCA plot,  PCA 3D plot.  
Dispersión del conjunto de datos:  plot de estimaciones de dispersión,  plot de media vs varianza.  
Resultados de expresión diferencial:  mapa de calor de expresión,  plot de expresión de top N genes,  plot de expresión emparejado de top N genes,  volcano plot,  volcano plot con etiquetas de genes significativos.  

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
</details>

**(*) <tablas_GSEA/ORA>.csv**  
GSEA para GO:BP, KEGG, Reactome;  
ORA de genes DOWN y UP para GO:BP, KEGG, Reactome
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
</details>

**(*) <figuras_GSEA/ORA>.png**  
GSEA para GO:BP, KEGG, Reactome: dotplot, ridgeplot, gseaplot2 de N top términos;   
ORA de genes DOWN y UP para GO:BP, KEGG, Reactome: dotplot, cnetplot(solo GO:BP), emapplot (solo GO:BP)
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
</details>


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
   | ANÁLISIS DE EXPRESIÖN DIFERENCIAL con DESeq2
   | (EXPLORACIÓN DE LOS DATOS)
   ├─→ [Figuras correlación entre muestras] ─→ heatmap correlacion, PCA, PCA 3D
   | (MODELAJE DE EXPRESIÓN DIFERENCIAL)
   ├─→ [Figuras dispersión de datos] ─→ estimaciones de dispersión, media vs varianza
   ↓ (CÁLCULO DE RESULTADOS)
A. [Métricas estadísticas para cada gen (baseMean, log2FC, lfcSE, stat Wald, pvalue, padj)]
B. [Lista de genes significativos UP y DOWN]
C. [Conteos normalizados]
   |   | VISUALIZACIÓN DE RESULTADOS
   .   ├─→ [Tabla de genes significativos con métricas estadísticas]
   .   └─→ [Figuras: Resultados de expresión diferencial] ─→ heatmap, volcano plot, plot de expresión, plot de expresión emparejado
   | ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL
   A ─→ |ranking stat| ─→ [GSEA: GO/KEGG/Reactome]
   B ─→ |UP/DOWN| ─→ [ORA: GO/KEGG/Reactome]
```


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
- Results_GSE132648_Souza2019: carpeta con todos los resultados sobre los datos de humano.  
- Results_GSE153648_Sorokin2023: carpeta con todos los resultados sobre los datos de ratón.  
- resultados_genes_significativos_SorokinSouza.xlsx : archivo Excel recopilatorio con todos los genes significativos del estudio.  
Dentro de ambas encontramos:  
    -  Carpeta "Enriquecimiento_Funcional": resultados de enriquecimiento (figuras y tablas; GSEA y ORA para GO:BP, KEGG, Reactome)  
    -  Figuras de correlación entre muestras (demuestran homogeneidad de las muestras):  
        -  heatmap de correlacion,  
        -  PCA,  
        -  PCA 3D.   
    -  Figuras de dispersión de datos (demuestran aptitud del modelo estadístico):  
        -  estimaciones de dispersión,  
        -  media vs varianza.    
    -  Figuras de resultados de expresión diferencial:  
        -  mapa de calor de expresión,  
        -  plot de expresión de top N genes,  
        -  plot de expresión emparejado de top N genes,  
        -  volcano plot,  
        -  volcano plot con etiquetas de genes significativos.    


## Scripts
Scripts de R.  
Estos scripts usan como entrada (input) los archivos de la carpeta Datasets y generan como salida (output) los archivos de la carpeta Results.  
- RNAseq_funciones.R : script con las funciones a utilizar.
- RNAseq_GSE132648_Souza2019.R : script para analizar datos de Souza et al. (humano).
- RNAseq_GSE153648_Sorokin2023.R : script para analizar datos de Sorokin et al. (ratón).


flowchart LR
[Datasets RNA-seq] --> [DESeq2] --> [Results] --> [Enriquecimiento]


flowchart LR
Counts[Conteos (GEO)] --> DESeq2
DESeq2 --> DEG[Genes DE]
DEG -->|UP/DOWN| ORA[ORA GO/KEGG/Reactome]
DESeq2 -->|ranking 'stat'| GSEA[GSEA GO/KEGG/Reactome]
GSEA --> Fig[Figuras y CSV]
ORA  --> Fig


