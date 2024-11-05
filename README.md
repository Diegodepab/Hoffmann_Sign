# Proyecto: signo de Hoffman

## To-Do List

- [X]  Ver si la hipótesis del trabajo es correcta o hay que mejorarla (esperar corrección)
- [X]  Objetivos (falta ver si hace falta un objetivo general y luego objetivos especificos, o si es una recopilación de objetivos)
- [X]  Materiales
    +  [X]  Human Phenotype Ontology
    +  [?]  R
    +  [X]  Python
          +  [ ] Librerías de python como iGraph, etc.
    +  [X] String
- [ ]  Método
    +  [X]  función que busca genes asociados a un fenotipo específico utilizando la API de **HPO**
          +  [ ]  Prepararlos en una URL los símbolos de los genes obtenidos para utilizarlos en la API de StringDB
    +  [ ]  Realizar una solicitud GET de la api **stringDB** (obteniendo una red de interacciones de proteínas)
          +  [ ] Decidir si filtrar por categorías o FDR (False Discovery Rate) 
    +  [X]  Evaluar la  la red de interacciones uso de **iGraph**
    +  [ ] Buscar más ideas como Buscar clusteres (grupos de genes conectados entre ellos), Usar Enriquecimiento Funcional para responder que hacen los genes en los clusteres? etc.
    proposiciones CHATGPT:
          +  [ ]  Algoritmos para hacer **clustering** ideales para este caso podrías explorar otros algoritmos como Louvain o Leiden para identificar comunidades o grupos de genes dentro de la red.
          +  [X] **Análisis Topológico de la Red**: Considera incluir métricas adicionales de evaluación de redes, como: Centralidad (Degree, Betweenness, Closeness). Coeficiente de agrupamiento (Clustering coefficient) y Densidad de la red
          +  [ ]  **Validación** del Modelo: Si tienes acceso a datos experimentales o publicaciones que puedan validar tus hallazgos, menciónalo en esta sección como una posible comparación o validación externa.
- [X] Entrega de la introducción (a espera de correcciones post evaluación)



## Preguntas para la tutoría/correo

- ¿El código debe ser un .sh obligatoriamente? 
- ¿Qué nivel de descarga debemos considerar, un equipo de ubuntu desde cero o simplemente instalar librerías necesarias?
- ¿Se podría usar un google Collab?

## Anotaciones generales

Este archivo servirá como referencia para el equipo. Aquí se anotará información útil para la búsqueda de datos y la organización general del proyecto.

- El **síndrome de Hoffman** no está relacionado con el **signo de Hoffman**.
- **Cuidado** con el **signo de Hoffman-Tinel**, ya que es diferente y puede generar confusión.
- Al buscar el **HPO** en **Monarch**, encontramos 12 enfermedades y 46 genes coincidentes con los de HPO.  
  - En **Orphanet**, no hay resultados relacionados.
  - En **OMIM**, aparecen 8000 entradas, pero no todas están vinculadas con el HPO.



