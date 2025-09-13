# TFM_bioinformatica

**Datasets**: Contiene los datos crudos sin editar, extraídos de GEO, que serán utilizados en los scripts.
- Souza et al. 2019: conteos de RNA-seq de muestras de sangre periférica humana, tras suplementación con aceite marino enriquecido en omega-3 o placebo.
- Sorokin et al. 2023: conteos de RNA-seq de muestras de ratón (aorta, hígado, piel), tras suplementación con DHA, EPA o placebo.

**Results**: Contiene los resultados obtenidos al aplicar los scripts.


**Scripts**: Contiene scripts de R; estos scripts usan como entrada los archivos de la carpeta Datasets y generan los archivos de la carpeta Results.
- RNAseq_funciones.R : script con las funciones a utilizar.
- RNAseq_GSE132648_Souza2019.R : script para analizar datos de Souza et al. (humano).
- RNAseq_GSE153648_Sorokin2023.R : script para analizar datos de Sorokin et al. (ratón).
