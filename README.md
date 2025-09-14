# TFM_bioinformatica

Este repositorio contiene conjuntos de datos de GEO procedentes de RNA-seq, scripts en R para análisis de expresión diferencial y enriquecimiento funcional, y los resultados generados.

Análisis de RNA-seq: preprocesado de conteos, DESeq2, y enriquecimiento funcional (GO/KEGG/Reactome).

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
[Datasets RNA-seq] --> [DESeq2] --> [Results] --> [Enriquecimiento]



<details>
Vista en árbol (resumida)


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
