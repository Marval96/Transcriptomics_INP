# **Entrenaminento en analisis de imagenes t-CyCIF**

Instructora: [Angela Zsabo](https://www.linkedin.com/in/ang%C3%A9la-szab%C3%B3-27321a125/)

Este analsis esta basado en el [repostorio](https://github.com/CruzOsuna/BMF_t-CyCIF/tree/main)  de Cruz y en el de [image_processing](https://github.com/CruzOsuna/BMF_t-CyCIF/tree/main) de [Färkkilä Lab](https://github.com/farkkilab). 



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

### **Registration: Stitching**

El codigo empleado para ejecutar este paso se basa en el repostorio de [farkkilab](https://github.com/farkkilab/image_processing/blob/main/pipeline/1_stitching/ashlar_workflow.py)
solo se modifcaron las rutas en nuestra computadora, muy similar a lo que indica el repositorio de Cruz.

        my_path = "/home/jrmarval/cycif/input/illimination_correction"
        output_path = "/home/jrmarval/cycif/output/registration"
        subfolders = [ f.path for f in os.scandir(my_path) if f.is_dir() ]
        file_type = 'czi' # 'nd2' 'rcpnl'

        illumination = 'Y' # 'N' 'Y'
        illumination_folder = "/home/jrmarval/cycif/output/illumination_correction"


Posteriormente ejecute el script ashlar_processing.py equivalente a stitching.py. Hice una variante para saber cuanto tarda en ejecutarse este proceso medinate un temporizador en un script de bash ejecutado en segundo plano:

        #!/bin/bash

        # Star timer
        start_time=$(date +%s.%N)

        # Run script ashlar
        python ashlar_processing.py -c 8

        # Execution time 
        end_time=$(date +%s.%N)
        execution_time=$(echo "$end_time - $start_time" | bc)

        # Convert seconds to minutes>seconds
        minutes=$(echo "scale=0; $execution_time / 60" | bc)
        seconds=$(echo "scale=0; $execution_time % 60" | bc)
        total_minutes=$(echo "$minutes + ($seconds > 0)" | bc)

        echo "Script Done"
        echo "Execution time: $minutes minutes $seconds seconds"
        echo

        # Alert by mail
        # Body mail
        BODY="Script done in: ${total_minutes} minutes."

        # Send mail
        echo "$BODY" | mail -s "Notice: Job Done" jhonatanraulm@gmail.com

        # End 
        echo "Done"

La estructura final del repositorio es la siguiente:

        ├── ashlar_processing.py
        ├── input
        │   └── illimination_correction
        │       └── BC577_20_1_5_cycles
        │           ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi
        │           ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi
        │           ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi
        │           ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi
        │           ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi
        │           └── BC577_20_1_041225_5_SLIDES_CICLO_6.czi
        ├── nohup.out
        ├── output
        │   ├── illumination_correction
        │   │   └── BC577_20_1_5_cycles
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi-dfp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_1.czi-ffp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi-dfp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_2.czi-ffp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi-dfp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_3.czi-ffp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi-dfp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_4.czi-ffp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi-dfp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_5.czi-ffp.tif
        │   │       ├── BC577_20_1_041225_5_SLIDES_CICLO_6.czi-dfp.tif
        │   │       └── BC577_20_1_041225_5_SLIDES_CICLO_6.czi-ffp.tif
        │   └── registration
        │       └── BC577_20_1_5_cycles.ome.tif
        └── run_ashlar.sh

El rsultado de este paso es el archivo: **BC577_20_1_5_cycles.ome.tif**