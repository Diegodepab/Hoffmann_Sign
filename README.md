# 🧬 **El Signo de Hoffmann como Indicador de Procesos Neurodegenerativos: Implicaciones Genéticas y Moleculares**
Este estudio analiza las interacciones entre genes y proteínas relacionadas con el **Signo de Hoffmann**, un indicador clínico asociado a enfermedades neurodegenerativas. Utilizando bases de datos bioinformáticas como la **Human Phenotype Ontology (HPO)** ([HP:0031993](https://hpo.jax.org/browse/term/HP:0031993)) y **StringDB**, se identificaron genes relevantes y se construyó una red de interacciones proteína-proteína. A partir de esta red, se aplicaron algoritmos de análisis de redes y enriquecimiento funcional para identificar grupos funcionales y vías biológicas clave. Los resultados revelaron la participación de procesos neuronales críticos, como la **plasticidad sináptica**, el **transporte intracelular** y la **regulación de la biosíntesis de ATP**. Este enfoque contribuye a una mejor comprensión de los mecanismos moleculares asociados al **Signo de Hoffmann**, ofreciendo posibles puntos de partida para futuras investigaciones en neurodegeneración.

## 🚀 **Resumen del Proyecto**

Este proyecto tiene como objetivo explorar y analizar interacciones genéticas relevantes para el **Signo de Hoffmann**, un fenotipo relacionado con enfermedades neurodegenerativas. A través de un conjunto de herramientas, se descargan datos de genes, se propagan interacciones mediante **DIAMOnD**, y se analizan propiedades de redes genéticas utilizando **R**.

Los pasos de este pipeline incluyen:

1. **Obtención de Datos Genéticos**: Se descargan datos de genes relacionados con el HPO.
2. **Conversión de Genes a IDs de STRING**: Los genes se convierten en identificadores compatibles con la base de datos **STRING**.
3. **Descarga y Preparación de la Red de STRING**: Se obtiene la red de interacciones de proteínas de STRING.
4. **Propagación de Genes con DIAMOnD**: Propagación de las interacciones a través de la red utilizando **DIAMOnD**.
5. **Análisis de la Red**: Se realizan análisis de la red con **R**.
6. **Enriquecimiento Funcional**: Se realiza un análisis de enriquecimiento funcional para identificar procesos biológicos relevantes.

---

## 🛠️ **Flujo de Trabajo (ejecutar `launch.sh`)

### El script `launch.sh` automatiza el siguiente proceso:

1. **Instalación de Dependencias**
   - Python y R son verificados e instalados si es necesario.
   - Paquetes de Python son instalados desde `requirements.txt`.
   - Paquetes de R son instalados desde `Rreqs.txt`.

2. **Descarga de Datos de Genes**
   - Se descargan datos de genes relacionados con el HPO desde una [fuente pública.](https://ontology.jax.org/api/network/annotation/HP:0031993/download/gene)

3. **Conversión de Genes a STRING IDs**
   - Los genes obtenidos son convertidos a identificadores STRING usando un script de Python (`genes2string.py`).

4. **Descarga de la Red de STRING**
   - La red de proteínas de STRING es [descargada](https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz) y descomprimida.

5. **Propagación de Genes con DIAMOnD**
   - Los genes mapeados son propagados a través de la red usando el algoritmo **DIAMOnD**.

6. **Obtención de Interacciones**
   - Interacciones adicionales de STRING son obtenidas mediante el script `string_interactions.py`.

7. **Análisis de la Red**
   - Se analizan las propiedades de la red usando R con el script `propiedades_red.R`.

8. **Análisis de Enriquecimiento Funcional**
   - Se realiza un análisis de enriquecimiento funcional para obtener insights sobre los procesos biológicos involucrados.

9. **Resultados**
   - Los resultados del análisis son generados y almacenados en la carpeta `results/`.

---

## 📦 **Requisitos**

Este pipeline requiere las siguientes dependencias (Podrás descargarlas mediante el `setup.sh`:

### 🐍 **Python**
- **NetworkX**
- **Requests**
- **Pandas**

### 📊 **R**
- **iGraph**
- **Bioconductor**


## 💡 **Cómo Usar Este Proyecto**

### 1. **Instalación de Dependencias**

Para instalar todas las dependencias necesarias, simplemente ejecuta el script `setup.sh`:

```bash
bash setup.sh
```
>[!info]
>Asegúrate de tener instalados Python y R en tu sistema para ejecutar correctamente el script `setup.sh.`

Esto instalará tanto las dependencias de Python como de R, configurando los entornos necesarios.

### 2. **Ejecutar el Pipeline**

```bash
bash launch.sh
```
Esto descargará los datos, ejecutará la propagación de genes, realizará el análisis de la red y generará los resultados en la carpeta results/.

## Autores

  _Autor:_ [martacuevasr](https://github.com/martacuevasr)

  _Autor:_ [Diegodepab](https://github.com/Diegodepab)

  _Autor:_ [JuanSoM](https://github.com/JuanSoM)
  
 _Autor:_ [AlexSilvaa9](https://github.com/AlexSilvaa9)
