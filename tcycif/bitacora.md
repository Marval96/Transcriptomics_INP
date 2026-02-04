# **Entrenaminento en analisis de imagenes t-CyCIF**

### Illumination Correction

Pasos previos: instalar docker y crear la imagen conforme a lo indicado en el 
[repositorio](https://github.com/CruzOsuna/BMF_t-CyCIF/tree/main/Processing_Stable/00_Illumination_correction) de Cruz. 

* Para evitar correr sudo en cada comando de Docker:

        sudo usermod -aG docker $USER

* Restaurar sesion para que se apliquen los cambios de Docker.

        newgrp docker

*  Verificar que Docker funciona

        docker run hello-world
    
* Construir la imagen de Docker **en la carpeta Scripts, aqui se encuentra el Docker file**

        docker build -t mybasic-image .

* Correr el contenedor modificando recursos computacionales:

        sudo docker run --privileged -it -m 20g --cpus=8 \
        --mount type=bind,source="$(pwd)",target=/data mybasic-image bash

* Dentro de contenedor ir a data:

    cd /data

Debe existir el directorio *input* y *output*. Las imagenes a analizar deberan estra en un subdirectorio en *input*. 

* Ejecutar el analsiis:

        bash BaSiC_run.sh

Al final el difrectorio con los resultados en *output* se ve asi:

    (base) jrmarval@BMFLAB:~/cycif/BMF_t-CyCIF/Processing_Stable/00_Illumination_correction$ tree
    .
    ├── IC_Verification
    │   ├── comparacion_marcadores_3_FA2664P53_2_5.png
    │   └── violyn_plot.py
    ├── Scripts
    │   ├── BaSiC_run.sh
    │   ├── Dockerfile
    │   ├── imagej_basic_ashlar.py
    │   ├── imagej_basic_ashlar_filepattern.py
    │   ├── input
    │   │   └── BC577_20_1_5_cycles
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi
    │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_6.czi
    │   │       └── output
    │   ├── live_progress_monitor.sh
    │   └── output
    │       └── BC577_20_1_5_cycles
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi-dfp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi-ffp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi-dfp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi-ffp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi-dfp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi-ffp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi-dfp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi-ffp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi-dfp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi-ffp.tif
    │           ├── BC577_20_1_041225_5_SLIDES_CICLO_6.czi-dfp.tif
    │           └── BC577_20_1_041225_5_SLIDES_CICLO_6.czi-ffp.tif
    ├── Seminar_image_processing_27022025.pdf
    ├── bitacora.md
    ├── file_format_example_1.png
    ├── file_format_example_2.png
    ├── readme.md
    └── run_commands