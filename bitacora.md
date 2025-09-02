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

---

### **Pendientes:**

+ Proporciones celulres !
+ Ensamblae transcriptómico de líneas celulares
+ Qué computadora utilizaré y cuáles son sus características? Si es Windows, puedo utlizar WSL con todos los recursos computacionales del equipo.
+ Parámetros experiemntales RNAseq Bulk

---

### **Literatura:**

+ [Spatial transcriptomics identifies molecular niche dysregulation associated with distal lung remodeling in pulmonary fibrosis](https://www.nature.com/articles/s41588-025-02080-x)

---




