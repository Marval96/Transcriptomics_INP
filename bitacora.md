# **Transcriptomics**

## Información general:

Bitácora de trabajo durante el periodo agosto 2025-enero 2026.

**Objetivo:** diseñar e implentar el flujos de trabajo bioinformáticos para el análsis de datos transcriptómicos. 

+ Contacto: Biol.Exp. Raúl Valderrama: jhonatanraulm@gmail.com

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

> Septimebre 5, 2025

**Métodos del artúiclo base para el análasis**





---

### RNAseq: FaDu

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

> Septiembre 4, 2025

**Análisis transcriptómico FaDU**

El trabajo del cual derivan estos datos de secuenciación es: [STING enhances cell death through regulation of reactive oxygen species and DNA damage](https://www.nature.com/articles/s41467-021-22572-8). Resumen: La **resistencia a los agentes que dañan el ADN es una causa significativa de fracaso del tratamiento** y malos resultados en oncología. Para identificar reguladores no reconocidos de supervivencia celular realizamos una pantalla CRISPR-Cas9 de genio entero usando el tratamiento con **radiación ionizante como presión selectiva**, e identificamos **STING (estimulador de genes de interferón) como un regulador intrínseco de supervivencia celular**. Mostramos que STING regula un programa transcripcional que **controla la generación de especies reactivas de oxígeno (ROS)**, y que la pérdida STING altera la homeostasis ROS para reducir el daño del ADN y causar resistencia terapéutica. De acuerdo con estos datos, el análisis de los tumores de los especímenes de pacientes con **carcinoma de células escamosas de cabeza y cuello** muestran que la expresión **STING baja se asocia con peores resultados**. También demostramos que la activación farmacológica de **STING mejora los efectos de la radiación ionizante** in vivo, proporcionando una justificación para las combinaciones terapéuticas de agonistas STING y agentes dañinos de ADN. Estos resultados destacan un papel de STING que está más allá de su función canónica en la **detección de dinucleotida cíclica y ADN**, e identifican a **STING como un regulador de la homeostasis celular ROS y la susceptibilidad de células tumorales a agentes reactivos dependientes del oxígeno.**

**Las células cancerígenas con poco STING generan menos radicales libres cuando se les daña el ADN, así que se dañan menos y resisten más la radioterapia. Si se activa STING con fármacos, las células tumorales acumulan más daño y la terapia funciona mejor.**

Análisis de expresión diferencial entre FaDu WT y Fadu irrWT:

Los punto de corte para el Fold Change es de 1.5 y p.adj < 0.05.

| Cambio en expresión (fold change) | Log2FC (valor a usar en GEO2R) | Interpretación |
|---------------------------------|-------------------------------|----------------|
| 1.5×                             | 0.58                          | Cambio moderado (50% más o menos) |
| 2×                               | 1                             | Duplicación o reducción a la mitad |
| 3×                               | 1.58                          | Triplicación o reducción a 1/3 |
| 4×                               | 2                             | Cuadruplicación o reducción a 1/4 |
| 8×                               | 3                             | Octuplicación o reducción a 1/8 |

El análisis diferencial identificó 366 genes diferencialmente expresados, de los cuales 82 estaban subexpresados y 284 sobreexpresados (Figura 1B). El análisis de enriquecimiento funcional de los genes sobreexpresados reveló asociación con procesos de estructura y remodelación de la matriz extracelular y señalización mediada por citocinas.

Posteriormente, se realizó un análisis de intersección entre tres conjuntos de datos: genes asociados a neo-loops, genes identificados por RNA-seq y genes correspondientes al listado de anticuerpos disponibles del grupo (Figura 3). Los resultados se resumen a continuación:

+ Neo-loops vs RNA-seq: tres genes compartidos, COL8A1, SPRR1A y SPRR1B.

+ Neo-loops vs Anticuerpos: un único gen compartido, COL1A1.

+ RNA-seq vs Anticuerpos: cuatro genes comunes TIMP3, FN1, MMP1, MMP2, y MMP13.

+ Intersección entre los tres conjuntos: No se identificaron genes compartidos

Todo el material relacionado se encuentra disponible [aquí](/rnaseq_fadu).

> Septiembre 5, 2025

Análisis de expresión diferencial entre **FaDu WT y Fadu IrrKo:**

Los punto de corte para el Fold Change es de 1.5 y p.adj < 0.05.

El análisis diferencial identificó 410 genes diferencialmente expresados, de los cuales 194 subexpresados y 216 sobreexpresados. El análisis de enriquecimiento funcional de los genes sobreexpresados reveló asociación con procesos de estructura y remodelación de la matriz extracelular, señalización mediada por citocinas y formación de la envoltura córnea.

Posteriormente, se realizó un análisis de intersección entre tres conjuntos de datos: genes asociados a neo-loops, genes identificados por RNA-seq y genes correspondientes al listado de anticuerpos disponibles del grupo (Figura 3). Los resultados se resumen a continuación:

+ Neo-loops vs RNA-seq: dos genes compartidos, SPRR1A y SPRR1B.

+ Neo-loops vs Anticuerpos: un único gen compartido, COL1A1.

+ RNA-seq vs Anticuerpos: tres genes comunes, MMP1, MMP2 y MMP13.

+ Intersección entre los tres conjuntos: No se identificaron genes compartidos

**Interpretación:** ambas comparativas de la línea FaDu, presentan una escasa intersección entre los conjuntos evaluados, lo que sugiere una superposición limitada de genes entre las diferentes fuentes de datos.

**Sugerencia:** dejar de lado el análisis con GEO2R... considero que podríamos obtener mejores resultados si manejamos los datos en crudo.

GEO2R es útil para **una exploración rápida** porque permite obtener genes diferencialmente expresados sin necesidad de programar. Tiene **limitaciones claras**:

* Trabaja sobre datos **ya procesados/normalizados**. 
* No ofrece control total sobre los **parámetros de normalización, filtrado o estadísticos**.

Sí nostros procesamos los datos algunas ventajas son:

* **Consistencia** en el preprocesamiento.
* **Flexibilidad** en los umbrales y métodos estadísticos.
* Posibilidad de incorporar **modelos más realistas**
* Resultados más **robustos y reproducibles**.


> + GEO2R = exploratorio, rápido, para validar ideas iniciales.
> + Datos crudos + análisis bioinformático = más tiempo y trabajo, pero resultados más confiables y con mayor aceptación académica.
> ¿Habrá diferencias? Sí, especialmente en el **número y tipo de genes identificados como diferencialmente expresados** y en la **resolución de los análisis de enriquecimiento**.


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




