# TFM_bioinformatica
Trabajo Fin de Máster del Máster en Bioinformática, Universidad Internacional de la Rioja (UNIR)
Autor: Víctor Saavedra Yturriagagoitia
Fecha: Septiembre 2025

Este repositorio contiene conjuntos de RNA-seq obtenidos de Gene Expression Omnibus (GEO), scripts en R para análisis de expresión diferencial (DESeq2) y enriquecimiento funcional (GO/KEGG/Reactome), y los resultados generados.



## Índice
- [Estructura del repositorio](#estructura-del-repositorio)
- [Requisitos](#requisitos)
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


<details>
<summary>Vista en árbol</summary>  

    
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
│  │  ├─ (*) <figuras_expresion_diferencial>.png
│  │  └─ Enriquecimiento_Funcional/
│  │     └─ Souza_aceitevsplacebo/
│  │        ├─ (*) <tablas_GSEA/ORA>.csv
│  │        └─ (*) <figuras_GSEA/ORA>.png
│  ├─ Results_GSE153648_Sorokin2023/
│  │  ├─ (*) resultados_genes_significativos__<tejido>_<trat1vstrat2>.csv
│  │  ├─ (*) <figuras_expresion_diferencial>.png
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
```
Nota: las ramas con (*) indican multiplicidad de archivos o carpetas.  

**(*) <figuras_expresion_diferencial>.png**  
    
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
```
* distinto para cada subset (tejido y combinación de tratamientos)  


**(*) <tablas_GSEA/ORA>.csv**  
(GSEA para GO:BP, KEGG, Reactome; ORA de genes DOWN y UP para GO:BP, KEGG, Reactome)  
     
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
```
* distinto para cada subset (tejido y combinación de tratamientos)  


**(*) <figuras_GSEA/ORA>.png**  
(GSEA para GO:BP, KEGG, Reactome: dotplot, ridgeplot, gseaplot2 de N top términos;   
 ORA de genes DOWN y UP para GO:BP, KEGG, Reactome: dotplot, cnetplot(solo GO:BP), emapplot (solo GO:BP))  
     
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
```
* distinto para cada subset (tejido y combinación de tratamientos)


</details>


---

## Flujo del análisis

```mermaid
flowchart LR
A[Conteos crudos (GEO)] --> B[DESeq2]
B --> C[Genes DE]
C -->|UP/DOWN| D[ORA (GO/KEGG/Reactome)]
B -->|ranking 'stat'| E[GSEA (GO/KEGG/Reactome)]
D --> F[Figuras y CSV en Results/]
E --> F
```


**Datasets**: Contiene los conteos crudos de RNA-seq sin editar, extraídos de GEO, que serán utilizados en los scripts; también los diseños experimentales.
- GSE132648_Souza2019: voluntarios humanos sanos suplementados con aceite marino enriquecido en omega-3 o placebo; muestras de sangre periférica.
    - archivo .txt: conteos crudos de RNA-seq
    - diseno....png: diseño experimental
- GSE153648_Sorokin2023: ratones suplementados con DHA, EPA o placebo; muestras de aorta, hígado, piel.
    - archivo .xlsx: conteos crudos de RNA-seq
    - diseno....png: diseño experimental

**Results**: Contiene los resultados obtenidos al aplicar los scripts.
- Results_GSE132648_Souza2019: carpeta con todos los resultados sobre los datos de humano.
    -  Enriquecimiento_Funcional: resultados de enriquecimiento funcional.
- Results_GSE153648_Sorokin2023: carpeta con todos los resultados sobre los datos de ratón.
    -  Enriquecimiento_Funcional: resultados de enriquecimiento funcional para cada pareja de tratamientos de cada tejido.
- resultados_genes_significativos_SorokinSouza.xlsx : archivo Excel con todos los genes significativos del estudio.

**Scripts**: Contiene scripts de R; estos scripts usan como entrada los archivos de la carpeta Datasets y generan los archivos de la carpeta Results.
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


