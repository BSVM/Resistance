# Installation and Usage of the Resistance Package

This document provides comprehensive instructions for installing and utilizing the `Resistance` R package for population genetic analysis of kdr mutations. The package standardizes statistical methods for calculating genotypic/allelic frequencies and implements classical population genetics analyses.

## Installation

### Step 1: Install `devtools` (if not present)
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

```
Step 2: Load devtoolss


```R
library(devtools)
```
Step 3: Install Resistance Package
mismo:


```R
devtools::install_github("BSVM/Resistance")

library(Resistance) 

```

## Practical Example: Analyzing kdr Mutation Dat

Below is a practical example demonstrating how to use key functions in the Resistance package.

# Example 1: Genotype Analysis of 12 Individuals
Crear un DataFrame con los Genotipos


This code creates a data frame where each row corresponds to an individual. The presence of a genotype is marked with "X" and absence with 0.

AA: Homozygous dominant

Aa: Heterozygous

aa: Homozygous recessive

```R

df <- data.frame(
  ID = 1:12,
  Pop = rep("Pop1", 12),
  AA = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
  Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
  aa = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

print(df)

```

Note: The package also provides a built-in dataset Loc1Pops3 with multiple populations, which can be accessed as follows:
```R

print (Loc1Pops3)

```
## Data Structure of Loc1Pops3s


| ID     | Pop   | AA   | Aa   | aa   |
|--------|-------|------|------|------|
| Cur 1  | Pop 1 | 0    | X    | 0    |
| Cur 2  | Pop 1 | 0    | X    | 0    |
| Cur 3  | Pop 1 | 0    | X    | 0    |
| Cur 4  | Pop 1 | 0    | X    | 0    |
| Cur 5  | Pop 1 | 0    | X    | 0    |
| Cur 6  | Pop 1 | 0    | X    | 0    |
| Cur 7  | Pop 1 | 0    | X    | 0    |
| Cur 8  | Pop 1 | 0    | X    | 0    |
| Cur 9  | Pop 1 | 0    | X    | 0    |
| Cur 10 | Pop 1 | 0    | X    | 0    |
| 11 filas   | ...   | ...  | ...  | ...  |


Para calcular las frecuencias Genotípicas (AA, Aa y aa) de la poblacion o poblaciones, utilizaran las funcion llamada `Import_single_loc`,
esta funcion hace un conteo de cada genotipo por poblacion, utilizando la informacion de presencia `X` y ausencia `O`.
Obtendremos una matriz llamda `df` que contiene el conteo de cuanto individuos presentan cada uno de los genotipos en la poblacion o poblaciones

```R
df1 <- Import_single_loc(Loc1Pops3)

print (df1)

```
La tabla `df1` tiene la siguiente estructura:

| Pop   | AA | Aa  | aa  |
|-------|----|-----|-----|
| Pop 1 | 4  | 33  | 3   |
| Pop 2 | 2  | 39  | 0   |
| Pop 3 | 7  | 22  | 1   |


A partir de la matriz de conteo `df1` , se calcularan las frecuencias genotípicas (AA, Aa y aa ) y alélicas (alelo A y alelo a), utilizando la función `Single_locus_frequencies`:

```R

df2 <- Single_loc_freq(df1)

print(df2)

```
La tabla `df2` tiene la siguiente estructura:
| Pop   |  n   |    aa      |    Aa     |    AA     |  Freq_A   |  Freq_a   |
|-------|------|------------|-----------|-----------|-----------|-----------|
| Pop 1 |  40  | 0.07500    | 0.82500   | 0.10000   | 0.51250   | 0.48750   |
| Pop 2 |  41  | 0.00000    | 0.95122   | 0.04878   | 0.52439   | 0.47561   |
| Pop 3 |  30  | 0.03333    | 0.73333   | 0.23333   | 0.60000   | 0.40000   |


Para realizar un test de Equilibrio de Hardy-Weinberg, utilizaremos la matriz de conteo de genotipos `df1` 

```R

resultados_HWQ <- HWQL1(df1)

print(resultados_HWQ)


```
Obtendran una tabla por todos los resultados, incluyendo el resultado del test chi quedrado con su respectivo p-valor

| Population | n  | aa     | Aa     | AA     | Freq_A | Freq_a | X2..df.       | P_Value |
|------------|----|--------|--------|--------|--------|--------|---------------|---------|
| Pop 1      | 40 | 0.07500 | 0.82500 | 0.10000 | 0.51250 | 0.48750 | 16.9537 (1)   | 0.00004 |
| Pop 2      | 41 | 0.00000 | 0.95122 | 0.04878 | 0.52439 | 0.47561 | 33.72688 (1)  | 0.00000 |
| Pop 3      | 30 | 0.03333 | 0.73333 | 0.23333 | 0.60000 | 0.40000 | 8.35648 (1)   | 0.00384 |


### Finalmente a partir de la matriz de datos de genotipos y alelos "df2" se puede graficar la frecuencia alelica por poblacion utilizando la funcion

- Primero deben reestructurar los datos de las frecuencias genotipicas y alelicas de las poblaciones en la matriz "df2" de manera que R pueda graficarlas
  para esto utilicen la funcion `RestructureAlleleFreq`. Esta funcion maneja dos argumentos de estrada, primero la matriz de datos y luego `Loc =` donde se debe
  indicar el numero de locis que tienen los datos. En este caso `Loc = 1` para un marcador monoloci. Si fuera un marcador bialelico es `Loc = 2`.


```R
df3 <- RestructureAlleleFreq(df2, Loc = 1)

```
 - Por ultimo podran graficar la frecuencia alelica se utilizara la siguiente funcion `PlotAlleleFrequency`, donde se debe indicar igualmente el numero de locis. `Loc = 1`
 
```R

 PlotAlleleFrequency(df3, Loc = 1)

```
<div align="center">
  <img src="https://github.com/user-attachments/assets/883f2b5a-f157-40ba-9b43-aee0d6574405" alt="Rplot" width="50%">
</div>


#

# Ejemplo de uso practico para analizar datos de mutacions KDR bialelicas

## Estructura de los datos originales

La tabla `Loc2Pops3` tiene la siguiente estructura:

```R

Print(Loc2Pops3)

```
## Frecuencias Genotípicas por Población

| ID      | Pop   | AA  | Aa  | aa  | BB  | Bb  | bb  |
|---------|-------|-----|-----|-----|-----|-----|-----|
| Cur 1   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 2   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 3   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 4   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 5   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 6   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 7   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 8   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 9   | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| Cur 10  | Pop 1 | 0   | X   | 0   | X   | 0   | 0   |
| ...     | ...   | ... | ... | ... | ... | ... | ... |


Para calcular las frecuencias Genotípicas (AABB, AABb, AAbb, AaBB, AaBb, Aabb, aaBB, aaBb y aabb) de la poblacion o poblaciones, utilizaran las funcion llamada `Import_di_loc`,
esta funcion hace un conteo de cada genotipo por poblacion, utilizando la informacion de presencia `X` y ausencia `O`.
Obtendremos una matriz llamda `dw1` que contiene el conteo de cuanto individuos presentan cada uno de los genotipos en la poblacion o poblaciones

```R

dw1 <- Import_di_loc(Loc2Pops3)

print(dw1)

```
La tabla `dw1` tiene la siguiente estructura:

## dw1$$MonGts

| Pop   | AA | Aa | aa | BB | Bb | bb |
|-------|----|----|----|----|----|----|
| Pop 1 |  4 | 33 |  3 |  3 |  1 |  1 |
| Pop 2 |  2 | 39 |  0 |  0 |  0 |  0 |
| Pop 3 |  7 | 22 |  1 |  1 |  7 |  0 |

## dw1$$DicGtps

| Pop   | AABB | AABb | AAbb | AaBB | AaBb | Aabb | aaBB | aaBb | aabb |
|-------|------|------|------|------|------|------|------|------|------|
| Pop 1 |    4 |    0 |    0 |   32 |    1 |    0 |    2 |    0 |    1 |
| Pop 2 |    2 |    0 |    0 |   39 |    0 |    0 |    0 |    0 |    0 |
| Pop 3 |    5 |    2 |    0 |   17 |    5 |    0 |    1 |    0 |    0 |


A partir de la matriz de conteo dw1$DicGtps, se calcularán las frecuencias alélicas (AB, Ab, aB y ab) utilizando la función `Diloc_Freq_Hd`. El cálculo haplotípico se realiza mediante un sistema de corrección que ajusta los haplotipos ambiguos generados por los genotipos heterocigotos AaBb.

El proceso comienza con una estimación inicial de las frecuencias haplotípicas, y cada iteración posterior refina estos cálculos utilizando una función de maximización. El número de iteraciones se define mediante el parámetro `max_iterations`, mientras que `tol = 1e-06` establece el nivel de tolerancia para la convergencia del modelo. Además, el argumento `dat` permite visualizar los resultados de las frecuencias haplotípicas, ya sea únicamente para la última iteración `dat = 1` o para todas las iteraciones en cada población `dat = 2`.

```R

dw2 <- Diloc_Freq_Hd(dw1$DicGtps, max_iterations = 10, tol = 1e-06, dat = 1)

prin(dw2)

```


### dw2 muestra el resultado de la frecuencia haplotipicas para cada poblacion, indicando en que iteracion el calculo converge y el valor de correcion utilizado.

| Iter     | Pop   | N  | AB       | Ab         | aB       | ab        | CorAB_ab   | CorAb_aB   |
|----------|-------|----|----------|------------|----------|-----------|------------|------------|
| 10       | Pop 1 | 40 | 0.5124999 | 1.048488e-07 | 0.4500001 | 0.03749990 | 0.9999975 | 2.454997e-06 |
| Initial  | Pop 2 | 41 | 0.5243902 | 0.000000e+00 | 0.4756098 | 0.00000000 | NA         | NA         |
| 10       | Pop 3 | 30 | 0.4992645 | 1.007355e-01 | 0.3840689 | 0.01593112 | 0.1705252 | 8.294748e-01 |

- Se debe reestructurar los datos de las frecuencias genotipicas y alelicas de las poblaciones en la matriz "dw2" de manera que R pueda graficarlas
  para esto utilicen la funcion `RestructureAlleleFreq`. Esta funcion maneja dos argumentos de estrada, primero la matriz de datos y luego `Loc =` donde se debe
  indicar el numero de locis que tienen los datos. En este caso `Loc = 2` para un marcador de dos loci. Si fuera un marcador monoalelico es `Loc = 1`.

```R

dw3 <- RestructureAlleleFreq(dw2, Loc = 2)

dw3

```

- Por ultimo podran graficar la frecuencia alelica se utilizara la siguiente funcion `PlotAlleleFrequency`, donde se debe indicar igualmente el numero de locis. `Loc = 2`
 
```R

 PlotAlleleFrequency(dw3, Loc = 2)

```
<div align="center">
  <img src="https://github.com/user-attachments/assets/01b2f176-e0a2-4ea7-a512-16345fa82768" alt="Rplot" width="50%">
</div>

