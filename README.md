# Instalación y uso del Paquete `Resistance`

Este documento proporciona instrucciones para instalar y utilizar el paquete `Resistance` desde GitHub. También incluye un ejemplo práctico de algunas de sus funciones principales.

## Instalación

### Paso 1: Instalar `devtools` (si no está instalado)

Asegúrese de tener instalado el paquete `devtools`. Puede verificarlo y, si es necesario, instalarlo ejecutando el siguiente código en R:

```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
```
### Paso 2: Cargar el paquete devtools
Cargue el paquete devtools para instalar paquetes desde GitHub:

```R
library(devtools)
```
### Paso 3: Instalar el paquete Resistance

Instale el paquete Resistance directamente desde GitHub y cargue el mismo:


```R
devtools::install_github("BSVM/Resistance") # Instalacion del paquete desde el depositorio github

library(Resistance) # Activacion del paquete

```

## Ejemplo de uso practico para analizar datos de mutacions KDR

A continuación, se presenta un ejemplo práctico que demuestra cómo utilizar algunas de las funciones del paquete Resistance
Para cargar datos de mutaciones kdr de un loci.

# Ejemplo 1: Análisis de Genotipos de 12 Individuos
Crear un DataFrame con los Genotipos

El siguiente código crea un data.frame que representa que genotipo presenta cada individuo con una "X", el "0" indica ausencia.
Los genotipos se representan como AA: Genotipo homocigoto dominante, Aa: Genotipo heterocigoto y aa: Genotipo homocigoto recesivo.
```R

df <- data.frame(
  ID = 1:12,
  Pop = rep("Pop1", 12),
  AA = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
  Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
  aa = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

print(df)

```

Nota: Opcionalmente el paquete proporciona una matriz de datos "Loc1Pops3" para multiples poblaciones, esta puede ser cargada utilizando: 
```R

print (Loc1Pops3) # Muestra la estructura de los datos originales

```

Calcular Frecuencias Genotípicas (AA, Aa y aa) de la poblacion o poblaciones, utilizaran las funcion llamada "Import_single_loc",
esta funcion hace un conteo de cada genotipo por poblacion, utilizando la informacion de presencia "X" y ausencia "O".
Obtendremos una matriz llamda "df" que contiene el conteo de cuanto individuos presentan cada uno de los genotipos en la poblacion o poblaciones

```R
df1 <- Import_single_loc(Loc1Pops3)

print (df1)

```
A partir de la matriz de conteo, se calcularan las frecuencias genotípicas (AA, Aa y aa ) y alélicas (alelo A y alelo a), utilizando la función Single_locus_frequencies:

```R

df2 < Single_loc_freq(df1)

print(df2)

```

| Pop   |  n   |    aa      |    Aa     |    AA     |  Freq_A   |  Freq_a   |
|-------|------|------------|-----------|-----------|-----------|-----------|
| Pop 1 |  40  | 0.07500    | 0.82500   | 0.10000   | 0.51250   | 0.48750   |
| Pop 2 |  41  | 0.00000    | 0.95122   | 0.04878   | 0.52439   | 0.47561   |
| Pop 3 |  30  | 0.03333    | 0.73333   | 0.23333   | 0.60000   | 0.40000   |


Para realizar un test de Equilibrio de Hardy-Weinberg, utilizaremos la matriz de conteo de genotipos df1.

```R
resultados_HWQ <- HWQL1(df1)

print(resultados_HWQ)


```
Obtendran una tabla por todos los resultados, incluyendo el resultado del test chi quedrado con su respectivo p-valor

### Finalmente a partir de la matriz de datos de genotipos y alelos "df2" se puede graficar la frecuencia alelica por poblacion utilizando la funcion

- Primero deben reestructurar los datos de las frecuencias genotipicas y alelicas de las poblaciones en la matriz "df2" de manera que R pueda graficarlas
  para esto utilicen la funcion "RestructureAlleleFreq". Esta funcion manera dos argumentos de estrada primero la matriz de datos y luego Loc = donde se debe
  el numero de loci que tienen los datos. En este caso Loc = 1. Si fuera un marcador bialelico, Loc = 2.


```R
df3 <- RestructureAlleleFreq(df2, Loc = 1)

```
 - Por ultimo podran graficar la frecuencia alelica utilizando la siguiente funcion
 
```R

 PlotAlleleFrequency(df3, Loc = 1)

```
<div align="center">
  <img src="https://github.com/user-attachments/assets/883f2b5a-f157-40ba-9b43-aee0d6574405" alt="Rplot" width="50%">
</div>

