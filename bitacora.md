# **Transcriptomics**

## Información general:

+ Contacto: Biol.Exp. Raúl Valderrama: jhonatanraulm@gmail.com

**Objetivo:** diseñar e implentar el flujos de trabajo bioinformáticos para el análsis de datos transcriptómicos. 

---

## Índice

- [General](#información-general)
- [ST: Lung](#spatial-transcriptomics-lung)
- [RNAseq bulk](#rnaseq-pablo)
- [Pendientes](#pendientes)
- [Literatura](#literatura)

---

### Spatial Transcriptomics: Lung 

> Septiembre 1, 2025

**Objetivo general: montar un flujo de trabajo para el análisis de datos transcriptómicos espaciales.**

+ Entragable 1: proporciones celulares, para mostrar en EUA.

Lectura y análisis del artículo [Spatial transcriptomics identifies molecular niche dysregulation associated with distal lung remodeling in pulmonary fibrosis](https://www.nature.com/articles/s41588-025-02080-x), para el trabajo de Miguel. 

Por hacer: 

+ Los datos estan libres? 
+ Qué plataforma se utiliza? XENIUM
+ Construir un set de datos de pruba, útiles para testear en flujo de trbaajo.
+ Cómo funciona Xenium? Cosnideraciones de la segmentación? Complejidad del clustering. 
+ Qué es un TMA?



---

### RNAseq: Pablo

> Septiembre 1, 2025

Para el análisis transcripómico, **lo ideal sería tener al menos triplicado biológico por comparación.** Sin embargo, se puede hacer con 2 replicas por condición, sacrificando el manejo de la variabilidad biológica por parte de DESeq2. 

> Al emplear una línea celular esta variabilidad se puede "compensar".

| Nivel de detección | Tipo de genes / objetivo                    | Réplicas por condición | Lecturas por muestra |
|-------------------|--------------------------------------------|-----------------------|--------------------|
| Bajo              | Cambios grandes y robustos                 | 2–3                   | 10–20 M            |
| Medio             | Genes medianamente expresados              | 3–4                   | 20–30 M            |
| Alto              | Genes poco expresados / isoformas          | 4–6                   | 40–50 M            |


**Notas:**

+ Más réplicas aumentan el poder estadístico, especialmente para genes de baja expresión

+ Para expresión diferencial “estándar” de genes codificantes, 20–30 M reads por muestra es suficiente.

+ Tip: se puede prorizar aumentar el número de replicas aunque se tengan menos reads.

> Septiembre 3, 2025

**¿Qué son las céulas FaDu?**

Es una línea celular humana derivada de un carcinoma de células escamosas de la hipofaringe. Las FaDu son conocidas por su crecimiento robusto y suelen emplearse en ensayos para comprender la proliferación de células cancerosas, la respuesta a agentes terapéuticos y la expresión de genes relacionados con la progresión del cáncer y la metástasis.

En la investigación científica, las células FaDu han sido fundamentales para examinar la eficacia de los tratamientos de radioterapia y quimioterapia, proporcionando información sobre las respuestas celulares al daño del ADN y los mecanismos de reparación. La versatilidad y relevancia de las células FaDu las convierten en un valioso modelo para la investigación oncológica, contribuyendo al desarrollo de terapias dirigidas y a la comprensión de la biología celular del cáncer a nivel molecular.

![FaDU](data_bitacora/FaDu.webp)

| Característica            | Valor       |
|----------------------------|-------------|
| Edad                       | 56 años     |
| Género                     | Hombre      |
| Etnia                      | Asiático    |
| Morfología                 | Epitelial   |
| Propiedades de crecimiento | Adherente   |

Las líneas celulares isogénicas son pares (o conjuntos) de líneas celulares genéticamente idénticas, excepto por una mutación específica introducida a propósito, lo que permite estudiar con claridad la función de ese gen o variante. 

En [GEO](https://www.ncbi.nlm.nih.gov/geo/) existe un transcriptoma RNAseq, de células FaDu isogénicas editadas con Cas9, en el cual analizan el efecto del KO de STING y radioterapia. El número de acceso es: [GSE147085](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147085). El conjunto de datos se localizo a través de [OmicsDI](https://www.omicsdi.org/).

Se hizo una análisis de expresión diferencial, con la plataforma GEO2R. El objetivo fue buscar en los DEG los genes de interés de la línea FaDu trabajada por el grupo, enfocados en genes asociados a la remodelación de matriz extracelular. Por los resultados encontrados, es muy probable que encontremos los genes de interés en el transcriptoma del grupo. 

**Nota:** Para tener una mejor simulación tal vez se podría hacer la coparativa entre el WT y el WT sometido a radiacíon.

+ ¿Cúal es el contexto experimental del grupo?¿Qué genes editaron?

**Tarea:** hacer la expresión diferencial entre WT y WT irradiado, de la línea FaDu. Es interesante saber que genes y vías se enriquecen tras recibir la lo dosis de radiacíon. La radioterapía genera daño genómico, lo mismo que sucede en AF, son dos condiciones similares que podrían orientarnos en el análisis de nuestro transcriptoma. 

+ Lo ideas sería tener datos FaDu+FA: buscar el set de datos
+ STING tiene relación con la misma clase de ruptura que induce FA? De ser así las comparativas y extrapolaciones serían más relevantes. 
+ El punto de corte de log2fc debe ser 1.5
+ Compartir resultados con Ulises y datos HI-C.
+ Explorar el enriquecimiento con ShinyGO.
---

### **Pendientes:**

+ Proporciones celulares !
+ Ensamblaje transcriptómico de líneas celulares
+ Qué computadora utilizaré y cuáles son sus características? Si es Windows, puedo utlizar WSL con todos los recursos computacionales del equipo.
+ Parámetros experimentales RNAseq Bulk
+ Cuál es el flujo de trabajo del NCBI en datos RNAseq bulk?

---

### **Literatura:**

+ [Spatial transcriptomics identifies molecular niche dysregulation associated with distal lung remodeling in pulmonary fibrosis](https://www.nature.com/articles/s41588-025-02080-x)

---




