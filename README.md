# TFM_bioinformatica

Análisis de RNA-seq para el TFM: preprocesado de conteos, DESeq2, y enriquecimiento funcional (GO/KEGG/Reactome).

> **Resumen**: Este repositorio contiene datasets desde GEO, scripts en R para análisis diferencial y enriquecimiento, y los resultados generados.

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
| **Datasets/** | Conjuntos de datos de RNA-seq de Gene Expression Omnibus (GEO). |
| **Results/**  | Figuras y tablas generadas por los scripts. |
| **Scripts/**  | Código en R creado en este trabajo. |


flowchart LR
A[Datasets] --> B[DESeq2]
B --> C[Results]
C --> D[Enriquecimiento]



<details>
<summary>Vista en árbol (resumida)</summary>


**Datasets**: Contiene los datos crudos sin editar, extraídos de GEO, que serán utilizados en los scripts.
- GSE132648_Souza2019: contiene conteos de RNA-seq de muestras de sangre periférica humana, tras suplementación con aceite marino enriquecido en omega-3 o placebo.
- GSE153648_Sorokin2023: contiene conteos de RNA-seq de muestras de ratón (aorta, hígado, piel), tras suplementación con DHA, EPA o placebo.

**Results**: Contiene los resultados obtenidos al aplicar los scripts.
- Results_GSE132648_Souza2019: carpeta con todos los resultados sobre los datos de humano.
- Results_GSE153648_Sorokin2023: carpeta con todos los resultados sobre los datos de ratón.
- resultados_genes_significativos_SorokinSouza.xlsx : archivo Excel con todos los genes significativos del estudio.

**Scripts**: Contiene scripts de R; estos scripts usan como entrada los archivos de la carpeta Datasets y generan los archivos de la carpeta Results.
- RNAseq_funciones.R : script con las funciones a utilizar.
- RNAseq_GSE132648_Souza2019.R : script para analizar datos de Souza et al. (humano).
- RNAseq_GSE153648_Sorokin2023.R : script para analizar datos de Sorokin et al. (ratón).
