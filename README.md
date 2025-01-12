# Instalación y Uso del Paquete `Resistance`

Este documento proporciona instrucciones para instalar y utilizar el paquete `Resistance` desde GitHub. También incluye un ejemplo práctico de algunas de sus funciones principales.

## Instalación

### Paso 1: Instalar `devtools` (si no está instalado)

Asegúrese de tener instalado el paquete `devtools`. Puede verificarlo y, si es necesario, instalarlo ejecutando el siguiente código en R:

```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
```
Paso 2: Cargar el paquete devtools
Cargue el paquete devtools para instalar paquetes desde GitHub:

R
Copiar código
library(devtools)
Paso 3: Instalar el paquete Resistance
Instale el paquete Resistance directamente desde GitHub con el siguiente comando:

R
Copiar código
devtools::install_github("usuario/Resistance")
Paso 4: Cargar el paquete Resistance
Después de instalarlo, cargue el paquete para comenzar a usarlo:

R
Copiar código
library(Resistance)
Ejemplo de Uso
A continuación, se presenta un ejemplo práctico que demuestra cómo utilizar algunas de las funciones del paquete Resistance.

Ejemplo 1: Análisis de Genotipos de 12 Individuos
Crear un DataFrame con los Genotipos
El siguiente código crea un data.frame que representa los genotipos de 12 individuos.

R
Copiar código
df <- data.frame(
  ID = 1:12,
  aa = c(0, 0, 0, 0, "X", 0, 0, 0, 0, 0, 0, "X"),
  Aa = c("X", 0, "X", "X", 0, "X", "X", "X", "X", "X", "X", 0),
  AA = c(0, "X", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)
Contar los Genotipos Presentes en la Muestra
Utilice la función Import_single_loc para generar una matriz de conteo basada en los datos:

R
Copiar código
df1 <- Import_single_loc(df)
Calcular Frecuencias Genotípicas y Alélicas
A partir de la matriz de conteo, calcule las frecuencias genotípicas y alélicas utilizando la función Single_locus_frequencies:

R
Copiar código
Single_locus_frequencies(df1)
Evaluar el Equilibrio de Hardy-Weinberg
Calcule el equilibrio de Hardy-Weinberg y realice una prueba de Chi-cuadrado para los datos:

R
Copiar código
resultados_HWQ <- calculate_HWQ(df1)
Imprimir los Resultados
Por último, imprima los resultados del análisis:

R
Copiar código
print(resultados_HWQ)
