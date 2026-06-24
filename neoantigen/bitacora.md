# **Deteccion de Neoantigenos**

La primera aproximacion a este analisis se hara con ENEO: [Efficient and effective identification of cancer neoantigens from tumor only RNA-seq](https://pmc.ncbi.nlm.nih.gov/articles/PMC12246784/). 

> ENEO permite identificar neoantígenos a partir de RNA-seq tumoral únicamente, utilizando un modelo bayesiano para compensar la ausencia de muestra normal, logrando un balance entre eficiencia, costo y precisión.

La identificación de neoantígenos tumorales es un paso clave en inmunoterapia personalizada, ya que estos péptidos derivados de mutaciones pueden ser reconocidos por el sistema inmune. El enfoque tradicional requiere múltiples tipos de datos:

+ Exoma tumor (WES)
+ Exoma normal (WES)
+ RNA-seq tumoral

Problemas:

+ Alto costo
+ Pipeline complejo
+ Baja eficiencia en algunos casos

ENEO, un método computacional que permite **detectar neoantígenos utilizando únicamente datos de RNA-seq tumoral**. El objetivo es reducir la complejidad y costo del pipeline sin perder precisión en la predicción. El principal desafío es diferenciar mutaciones somáticas (tumorales) de variantes germinales (tejido sana) sin una muestra nde de tejido normal.

**ENEO propone un modelo bayesiano basado en inferencia probabilística para estimar la probabilidad de que una variante sea:**

+ Somática (tumor)
+ Germinal (normal)
+ Error técnico

Utilizada datos de RNA-seq y bases de datos de variantes conocidas asi como patrones de expresión.

El pipeline de ENEO en terminos genrales es asi: 

1. Llamado de variantes (variant calling)
Identificación de SNPs e indels directamente desde RNA-seq
2. Filtrado probabilístico
Eliminación de ruido técnico
Clasificación de variantes mediante modelo bayesiano
3. Generación de péptidos mutados
Traducción de variantes a secuencias proteicas
4. Predicción de unión a MHC
Evaluación de afinidad de péptidos a moléculas HLA.

El resultado es una lista de neoantígenos candidatos.

El trabajo se valido em un dataset principal, TESLA benchmark (referencia en estudios de neoantígenos). Una validación adicional de 
2 cohortes independientes sobre diferentes tipos tumorales y distintos protocolos experimentales.

Ventajas: 

+ la identificación eficiente de neoantígenos derivados de mutaciones en DNA es comparables a métodos estándar.
+ Menor costo
+ Menor tiempo de análisis
+ Menor requerimiento de datos
+ Permite detectar: mutaciones específicas de RNA, RNA editing y Splicing aberrante. Esto amplía el repertorio potencial de neoantígenos

Limitaciones (No reemplaza completamente WES, pero es una alternativa eficiente):

+ Dependencia de la expresión génica;mutaciones en genes poco expresados pueden no detectarse.
+ Mayor ruido en variant calling desde RNA-seq
+ Cobertura desigual del transcriptoma

Este enfoque demuestra que:

+ Es posible realizar predicción de neoantígenos de forma más accesible utilizando únicamente RNA-seq, manteniendo precisión razonable.

**Notas para implementación:**
La cobertura (read depth) en RNA-seq es crítica para detectar variantes
Priorizar genes altamente expresados,
control de calidad estricto en variant calling, integrar con herramientas de predicción de MHC.

La identificación de neoantígenos tradicional depende de datos de DNA (WES/WGS) y RNA-seq, pero este enfoque presenta limitaciones importantes, incluyendo la sobreestimación de candidatos y baja correlación con péptidos realmente presentados por MHC.

Evidencia reciente indica que el transcriptoma (RNA-seq) representa mejor el repertorio de neoantígenos, al capturar no solo mutaciones expresadas, sino también alteraciones específicas de RNA como splicing aberrante y RNA editing.

Sin embargo, el uso exclusivo de RNA-seq enfrenta desafíos técnicos, como la ausencia de un control normal y un mayor nivel de ruido en el llamado de variantes.

ENEO aborda estas limitaciones mediante un modelo bayesiano que permite inferir la naturaleza somática de las variantes usando únicamente datos tumorales, habilitando así un pipeline completo de predicción de neoantígenos basado exclusivamente en RNA-seq.

Notas:  Mutaciones no sinónimas: cambios en el DNA que alteran la secuencia de aminoácidos de una proteína, permitiendo la generación de neoantígenos.

- MHC (Major Histocompatibility Complex): conjunto de moléculas encargadas de presentar péptidos al sistema inmune.

- HLA (Human Leukocyte Antigen): versión humana del MHC.

- MHC clase I:
  - Presenta péptidos intracelulares
  - Activa linfocitos T CD8+
  - Principal vía en reconocimiento de células tumorales

- MHC clase II:
  - Presenta péptidos extracelulares
  - Activa linfocitos T CD4+
  - Asociado a respuesta inmune adaptativa auxiliar

### **Flujo de trabajo**

+ Crear ambinete conda para la instlacion y gestion de la herramienta

    conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake apptainer

+ Clonar repositorio

    git clone https://github.com/ctglab/ENEO.git

**Configurar setup**
ENEO heavily depends on public genetic databases for germline probability estimation plus other fairly common resources daily needed for bioinformatics. In order to make the download and configuration of resources for the workflow as smooth as possible, **a python configuration script is available inside setup/download_res.py**

+ Crear ambiente conda para gestionar la configuracion:

    conda env create -f setup_env.yml

El archivo original causaba muchos conflictos con Pip, la version de Python, depdendcias etc...entonces se hizo una version simplificada:

    name: eneo_setup
    channels:
    - conda-forge
    - bioconda
    - defaults
    dependencies:
    - python=3.10
    - cyvcf2=0.31.1
    - bcftools
    - bedtools
    - samtools
    - gatk4
    - pip
    - pip:
        - certifi==2024.8.30
        - charset-normalizer==3.4.0
        - click==8.1.7
        - coloredlogs==15.0.1
        - humanfriendly==10.0
        - idna==3.10
        - numpy==2.1.2
        - pandas==2.2.3
        - python-dateutil==2.9.0.post0
        - pytz==2024.2
        - pyyaml==6.0.2
        - requests==2.32.3
        - rich==14.0.0
        - six==1.16.0
        - tzdata==2024.2
        - urllib3==2.2.3

Los archivos necesarios estan en: /home/jrmarval/neoantigen/ENEO/resources.

El archivo de configuracion quedo asi:

    OUTPUT_FOLDER: /home/jrmarval/neoantigen/ENEO/eneo_output/
    TEMP_DIR: /home/jrmarval/neoantigen/ENEO/eneo_temp/
    datadirs:
    BQSR: BQSR
    HLA_typing: HLA_typing
    VCF: VCF
    VCF_out: VCF_out
    bams: bams
    expression: expression_data
    index_folder: genome_index
    logs:
        align: log/align
        annotate_variants: log/annotate_variants
        bam_cleaning: log/bam_cleaning
        bam_readcount: log/bam_readcount
        base_recalibration: log/base_recalibration
        decompose: log/decompose
        export_quant: log/export_quant
        intervals: log/intervals
        pMHC: log/pMHC
        salmon_quant: log/salmon_quant
        snv_calling: log/snv_calling
        sort_bam: log/sort_bam
        star_idx: log/star_idx
        t1k: log/t1k
        trimming: log/trimming
    mapped_reads: mapped_reads
    peptides: peptides
    salmon_idx: salmon_index
    salmon_quant: quantification
    trimmed_reads: trimmed_reads
    trimming_report: fastp_report
    utils: utils
    execution_mode: full
    params:
    BQSR:
        RAM: 30000
        threads: 4
    deepvariant:
        threads: 4
        extra: "split_skip_reads=true,channels=''"
    gatk:
        RAM: 20
        extra:
        RGPU: unit1
        RGSM: 20
    MarkDuplicates:
        RAM: 30000
        threads: 4
    pMHC:
        calibration_frame: workflow/supplementary_res/optimal_percentile_netmhcpan.csv
        filter_peptides_script: workflow/scripts/filter_peptides.py
        germProb: 0.5
        hla_ligand_atlas: workflow/supplementary_res/HLA_ligand_atlas.tsv.gz
        max_length: 12
        min_length: 8
        netmhcpan_launcher_script: workflow/scripts/netmhcpan_launcher.py
        threads: 4
    STAR:
        RAM: null
        extra: '--twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat '
        threads: 12
    SplitNCigarReads:
        RAM: 30000
        threads: 4
    salmon:
        RAM: null
        extra:
        extra: --gcBias --seqBias --reduceGCMemory
        index: --keep-duplicates
        libtype: A
        zip_ext: gz
        threads: 8
    samtools:
        threads: 4
    strelka2:
        threads: 8
    t1k:
        dat_file: workflow/supplementary_res/hla.dat
        threads: 8
    vcfanno:
        threads: 8
        toml_script: workflow/scripts/createTOML.py
        vcfanno_binary: workflow/utils/vcfanno_linux64
        vcfanno_lua: workflow/utils/custom.lua
        vcfanno_toml: workflow/utils/vcfanno.toml
    vep:
        extra:
        assembly: GRCh38
        filtering: --gencode_basic --pick --af --check_existing --coding_only --format
            vcf --vcf --symbol --terms SO --no_intergenic --tsl
        plugins:
            Frameshift: workflow/utils/vep_plugins/Frameshift.pm
            Wildtype: workflow/utils/vep_plugins/Wildtype.pm
    # resources:
    #   dbsnps: path/to/dbsnps_withAF.vcf.gz
    #   deepvariant_rna_model: path/to/deepvariant_rna_model
    #   genome: path/to/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
    #   germline_prob_script: workflow/scripts/germProb.py
    #   giab_intervals: workflow/supplementary_res/GRCh38_giab_merged.bed.gz
    #   gnomad: path/to/af-only-gnomad.hg38.vcf.gz
    #   gsnps: path/to/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    #   gtf: path/to/gencode.v47.primary_assembly.annotation.gtf
    #   hla_script: workflow/scripts/HLA_typing.py
    #   indel: path/to/Homo_sapiens_assembly38.known_indels.vcf.gz
    #   intervals_coding: workflow/supplementary_res/intervals_coding.BED.gz
    #   PoN: path/to1000g_pon.hg38.vcf.gz
    #   REDI: path/to/REDI.BED.gz
    #   rna_errors_script: workflow/scripts/filter_rna_errors.py
    #   t1k_file: workflow/supplementary_res/hlaidx_rna_seq.fa
    #   toml_script: workflow/scripts/createTOML.py
    #   transcriptome: path/to/gencode.v47.transcripts.fa
    #   vep_cache: path/to/homo_sapiens_vep_105_GRCh38
    # slurm_log_dir: slurm-logs

    resources:
    dbsnps: /home/jrmarval/neoantigen/ENEO/resources/dbsnps_withAF.vcf.gz
    deepvariant_rna_model: /home/jrmarval/neoantigen/ENEO/resources/deepvariant_rna_model
    genome: /home/jrmarval/neoantigen/ENEO/resources/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
    germline_prob_script: workflow/scripts/germProb.py
    giab_intervals: workflow/supplementary_res/GRCh38_giab_merged.bed.gz
    gnomad: /home/jrmarval/neoantigen/ENEO/resources/af-only-gnomad.hg38.vcf.gz
    gsnps: /home/jrmarval/neoantigen/ENEO/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    gtf: /home/jrmarval/neoantigen/ENEO/resources/gencode.v45.annotation.gtf.gz
    hla_script: workflow/scripts/HLA_typing.py
    indel: /home/jrmarval/neoantigen/ENEO/resources/Homo_sapiens_assembly38.known_indels.vcf.gz
    intervals_coding: workflow/supplementary_res/intervals_coding.BED.gz
    PoN: null
    REDI: null
    rna_errors_script: workflow/scripts/filter_rna_errors.py
    t1k_file: workflow/supplementary_res/hlaidx_rna_seq.fa
    toml_script: workflow/scripts/createTOML.py
    transcriptome: /home/jrmarval/neoantigen/ENEO/resources/Homo_sapiens.GRCh38.cdna.all.fa.gz
    vep_cache: /home/jrmarval/neoantigen/ENEO/resources/vep_cache

**Ejecucion:**
    snakemake \
    --use-singularity \
    --cores 8 \
    --singularity-args "-B /home/jrmarval/neoantigen:/home/jrmarval/neoantigen"

Si alogo falla, ENEO puede volver a correr desde el punto de error gracias a su desarrollo modular:

    snakemake --use-singularity --cores 8  --bind /home/jrmarval/neoantigen:/home/jrmarval/neoantigen --rerun-incomplete
---
    snakemake --cores  8--use-singularity \
--singularity-args "--bind /home/jrmarval/neoantigen:/home/jrmarval/neoantigen"

**Comando funcinal con base en CODEX:**

    cd /home/jrmarval/neoantigen/ENEO
---

    env XDG_CACHE_HOME=/home/jrmarval/neoantigen/ENEO/.cache TMPDIR=/home/jrmarval/neoantigen/ENEO/eneo_temp conda run -n snakemake snakemake -s workflow/Snakefile --use-apptainer --apptainer-args '--unsquash -B /home/jrmarval/neoantigen:/home/jrmarval/neoantigen' --cores 1 --printshellcmds --rerun-incomplete --rerun-triggers mtime

## **Analsis de Expresion Diferencial**

#### **Control de calidad de datos RNA-seq**

Antes de iniciar el análisis de expresión diferencial se realizó control de calidad sobre los mismos archivos FASTQ usados para la detección de neoantígenos. El objetivo de este paso es evaluar la calidad de las lecturas crudas de RNA-seq, identificar posibles problemas técnicos y generar un reporte integrado que facilite la revisión de todas las muestras antes de continuar con el alineamiento, ensamblado, cuantificación y análisis de expresion diferencial.

El análisis se hizo con archivos paired-end (`R1` y `R2`) ubicados en:

```bash
/home/jrmarval/neoantigen/fastq_files/raw_data
```

El ambiente conda utilizado fue:

```bash
conda activate QualityControl
```

Versiones registradas desde el ambiente `QualityControl`:

| Herramienta | Versión |
|---|---:|
| FastQC | v0.12.1 |
| MultiQC | 1.31 |

Archivos FASTQ analizados:

| Muestra | Lectura R1 | Tamaño R1 | Lectura R2 | Tamaño R2 | Tipo de datos |
|---|---|---:|---|---:|---|
| FAHNSCC5 | `FAHNSCC5_S1_R1_001.fastq.gz` | 7.9G | `FAHNSCC5_S1_R2_001.fastq.gz` | 8.0G | RNA-seq paired-end |
| FASCE12 | `FASCE12_S2_R1_001.fastq.gz` | 11G | `FASCE12_S2_R2_001.fastq.gz` | 12G | RNA-seq paired-end |
| IC242054_1 | `IC242054_1_S3_R1_001.fastq.gz` | 2.8G | `IC242054_1_S3_R2_001.fastq.gz` | 3.0G | RNA-seq paired-end |
| IC242054_2 | `IC242054_2_S4_R1_001.fastq.gz` | 1.9G | `IC242054_2_S4_R2_001.fastq.gz` | 2.0G | RNA-seq paired-end |
| IC242054_3 | `IC242054_3_S5_R1_001.fastq.gz` | 1.2G | `IC242054_3_S5_R2_001.fastq.gz` | 1.3G | RNA-seq paired-end |

El control de calidad se ejecutó con dos scripts:

- `interaccion.sh`: script interactivo que verifica que el ambiente `QualityControl` exista y que `fastqc` y `multiqc` estén disponibles. Si el usuario confirma la ejecución, lanza `analisis.sh` en segundo plano con `nohup`.
- `analisis.sh`: script principal que ejecuta FastQC sobre cada archivo `*.f*q.gz`, guarda los reportes individuales en `FastQC_Results/` y después ejecuta MultiQC para generar un reporte global.

Código usado para iniciar el análisis:

```bash
#!/bin/bash

echo
echo "Hola, $USER !"
echo
echo -e "Este script realiza un análisis de calidad de datos obtenidos mediante secuenciación de nueva generación (NGS) de tipo paired-end.\n\
Utiliza la herramienta FastQC y MultiQC.\n\
Primero realiza el análisis de calidad para cada una de las muestras con FastQC y después con MultiQC genera un resumen global más fácil de interpretar."
echo
echo -e "Recuerda activar el ambiente Conda llamado 'QualityControl' el cual contiene las herramientas necesarias para el análisis.\n\
Para tener el ambiente Conda, sigue los siguientes pasos:\n\
    1. Descargar el archivo QualityControl.yml, el cual contiene el ambiente Conda.\n\
    2. Generar el ambiente Conda ejecutando: conda env create -f QualityControl.yml\n\
    3. Activar el ambiente Conda ejecutando: conda activate QualityControl."
echo
echo "¿Estás listo para ejecutar el análisis? (si/no)"
read respuesta
echo

if [[ "$respuesta" == "si" || "$respuesta" == "SI" || "$respuesta" == "Si" ]]; then
    if conda info --envs | grep -q "QualityControl" && command -v fastqc >/dev/null && command -v multiqc >/dev/null; then
        echo "¡Perfecto! El ambiente 'QualityControl' está activo y contiene FastQC y MultiQC."

        echo "Comenzando el análisis..."
        nohup ./analisis.sh &
        pid=$!
        echo "El análisis ha comenzado en segundo plano. El PID es: $pid"
        echo "El registro se guardará en el archivo 'nohup.out'. Puedes monitorearlo ejecutando: tail -f nohup.out"
    else
        echo "Error: El ambiente 'QualityControl' no está activo o no contiene FastQC y/o MultiQC."
        echo "Por favor, asegúrate de que el ambiente esté activo y que las herramientas estén instaladas."
        exit 1
    fi
elif [[ "$respuesta" == "no" || "$respuesta" == "NO" || "$respuesta" == "No" ]]; then
    echo "Gracias, vuelva pronto, mil besos."
    exit 0
else
    echo "Respuesta no válida. Por favor, responde con 'si' o 'no'."
    exit 1
fi
echo
```

Código usado para ejecutar FastQC y MultiQC:

```bash
#!/bin/bash

start_time=$(date +%s.%N)

echo "Iniciando análisis de calidad con FastQC para cada archivo..."
mkdir -p FastQC_Results

for file in *.f*q.gz; do
    echo "Procesando $file ..."
    fastqc -o FastQC_Results "$file"
done

echo
echo "Generando informe global con MultiQC..."
multiqc FastQC_Results -o FastQC_Results

end_time=$(date +%s.%N)
execution_time=$(echo "$end_time - $start_time" | bc)

minutes=$(echo "scale=0; $execution_time / 60" | bc)
seconds=$(echo "scale=0; $execution_time % 60" | bc)
total_minutes=$(echo "$minutes + ($seconds > 0)" | bc)

echo "Análisis de calidad completado. Los resultados se encuentran en la carpeta FastQC_Results."
echo "Tiempo de ejecución: $minutes minutos $seconds segundos."
echo

BODY="El análisis de calidad se hizo en ${total_minutes}."

echo "$BODY" | mail -s "Aviso: FA-HNSCC Control de calidad" jhonatanraulm@gmail.com

echo "Adjunto el reporte completo de MultiQC." | \
    mail -A "FastQC_Results/multiqc_report.html" -s "Reporte MultiQC: Control de calidad" jhonatanraulm@gmail.com

echo "Done"
```

Comandos útiles para reproducibilidad y monitoreo:

```bash
cd /home/jrmarval/neoantigen/fastq_files/raw_data
conda activate QualityControl
bash interaccion.sh
tail -f nohup.out
```

Los resultados esperados del control de calidad son:

- Reportes individuales de FastQC para cada archivo FASTQ en `FastQC_Results/`.
- Reporte global de MultiQC en `FastQC_Results/multiqc_report.html`.
- Registro de ejecución en `nohup.out` cuando el análisis se lanza en segundo plano.

