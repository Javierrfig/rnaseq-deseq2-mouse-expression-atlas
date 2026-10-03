# Análisis de expresión diferencial con DESeq2 (Expression Atlas E-MTAB-10413, *Mus musculus*)

Este repositorio contiene mi primer flujo de análisis de RNA-seq para datos de Expression Atlas del experimento E-MTAB-10413, que compara ratones knockout para *Slc38a10* con ratones wild-type en 4 condiciones de estrés ambiental y tiempo:
1. Sin tratamiento, tiempo 0 horas (none at 0 hour).
2. Privación de B27, 8 horas (B27 starved at 8 hour).
3. Reposición de aminoácidos, 1 hora (amino acid refed at 1 hour).
4. Privación de aminoácidos, 2 horas (amino acid starved at 2 hour).

## Contexto

Este proyecto replica y extiende un tutorial de análisis de expresión diferencial con DESeq2. Los datos originales provienen de [Expression Atlas](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10413/Results).

## Métodos

1. Descarga de datos de conteos y metadatos desde Expression Atlas.
2. Preprocesamiento de datos (filtrado, normalización, visualización de muestras).
3. Análisis de expresión diferencial con DESeq2:
   - Análisis A: Modelo con todas las condiciones juntas (`~genotype`), para una visión global.
   - Análisis B: Modelo por separado para cada condición, para replicar los resultados de Expression Atlas.
4. Validación de resultados del análisis B con genes de interés reportados en Expression Atlas.
5. Análisis de enriquecimiento de genes diferencialmente expresados en la vecindad del gen knockout *Slc38a10*, usando la prueba exacta de Fisher.
6. Análisis de sensibilidad para el umbral de log2 fold change usado para definir genes "alterados" en el análisis de enriquecimiento.

Detalles completos en el script `deseq2_differential_expression.R`.

## Resultados

- El gen knockout (Slc38a10) muestra una caída de expresión consistente en las 4 condiciones (log2FC de -6.3 a -7.5), confirmando que el knockout funcionó correctamente (ver boxplot de Slc38a10 generado en RStudio).
- Los resultados del Análisis B concuerdan con los de Expression Atlas. Por ejemplo, Arsi en la condición "none at 0 hour" da un log2FC de -9.14 en R y -9.1 en Expression Atlas, p ajustado: 4.43e-17 vs 3.55e-14 (ver deseq_por_condicion.csv).
- Los genes que pasan el filtro estricto de Expression Atlas (|log2FC| ≥ 7 y p ajustado ≤ 0.06) coinciden con los reportados en la interfaz web (ver tabla_tipo_gxa.csv), incluyendo la distribución por condición: Col8a1, Madcam1, Anxa1 y A730049H05Rik aparecen en "B27 starved at 8 hour", mientras que Arsi, Glycam1 y Acta2 aparecen en "none at 0 hour".
- No se encontró evidencia de enriquecimiento de genes alterados en la vecindad de Slc38a10 (53-54 genes en la ventana de ±1 Mb) en ninguna de las 4 condiciones, p ajustado ≥ 0.5, prueba exacta de Fisher con |log2FC| ≥ 1 (ver vecinos_gen_knockout.csv y figures/cromosoma_gen_knockout.png). El análisis de sensibilidad confirmó que esta conclusión se mantiene con los 4 umbrales probados (ver vecinos_sensibilidad.csv y figures/sensibilidad_vecinos.png).

## Discusión

Este flujo de trabajo replica exitosamente un análisis estándar de RNA-seq utilizando DESeq2 [2]. La concordancia de los resultados con Expression Atlas [1] valida la implementación del pipeline.

La principal diferencia respecto al tutorial original es el diseño experimental. En E-MTAB-5244, una sola comparación con ~genotype basta para replicar Expression Atlas. En E-MTAB-10413, el factor condición (estrés ambiental × tiempo) introduce variabilidad que el modelo simple no distingue: al mezclar muestras de "none at 0 hour" con muestras de "B27 starved at 8 hour", DESeq2 estima un efecto promedio que no corresponde a ninguna comparación de Expression Atlas. Separar el análisis por condición resuelve esto. Comprender cómo el diseño del experimento determina la fórmula del modelo es un concepto clave en el análisis de expresión diferencial [3].

El análisis exploratorio de genes vecinos al locus de Slc38a10 no encontró evidencia de un efecto local. Estudios previos han documentado que en ratones knockout generados con células madre embrionarias de la cepa 129, los genes vecinos al locus editado pueden mostrar expresión diferencial debido a mutaciones pasajeras heredadas de la cepa donante, incluso tras más de 10 generaciones de retrocruce [4][5]. En este caso, la ausencia de enriquecimiento podría deberse a diferencias en el método de generación del knockout o en el contexto biológico. Un análisis de sensibilidad con 4 umbrales distintos de log2 fold change confirmó que esta conclusión no depende de la elección arbitraria de un umbral.

Como limitaciones, este análisis no exploró otros tamaños de ventana genómica ni otros métodos de detección de efectos locales. Además, con 4 réplicas por grupo y condición, la potencia estadística es limitada para detectar efectos sutiles.

## Referencias

[1] Papatheodorou I, et al. (2020). "Expression Atlas update: from tissues to single cells." Nucleic Acids Research, 48(D1), D77-D83. doi: 10.1093/nar/gkz947.

[2] Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

[3] Conesa A, et al. (2016). "A survey of best practices for RNA-seq data analysis." Genome Biology, 17, 13. doi: 10.1186/s13059-016-0881-8.

[4] Vanden Berghe T, et al. (2015). "Passenger Mutations Confound Interpretation of All Genetically Modified Congenic Mice." Immunity, 43(1), 200-209. doi: 10.1016/j.immuni.2015.06.011.

[5] Westphal D, et al. (2017). "Passenger gene confounds phenotypes of SARM1-deficient mice." (preprint: bioRxiv doi: 10.1101/2021.08.25.457655).

Enlaces adicionales:
- Análisis de la expresión diferencial de ARN-Seq: un proyecto de bioinformática - OMGenomics
- DESeq2 Design Explained: Understanding Interaction Terms in RNA-Seq Analysis - Bioinformagician
- Expression Atlas E-MTAB-10413 Results
