#' Funcion para calcular la frecuencias haplotipica con mutaciones Dilocus
#'
#' @param data Datos originales se puede utilizar V1016I.F1535C o V1016G.F1535C
#' @param pais Pais del cual deseamos extraer los datos
#'
#' @return Tabla con las frecuencias haplotipicas del pais por año
#'
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#'
#' @examples
#' Freq_hap_mexico <- Hap_DilocusfreqforYear(V1016I.F1534C, "Mexico")
#'
#' @export
Hap_DilocusfreqforYear <- function(data, pais) {
  datos_filtrados <- data %>%
    filter(Country == pais)

  # Calcular las frecuencias promedio por año
  frecuencias_promedio <- datos_filtrados %>%
    group_by(`Year of Sample`) %>%
    summarise(
      R3 = mean(R3, na.rm = TRUE),
      R2 = mean(R2, na.rm = TRUE),
      R1 = mean(R1, na.rm = TRUE),
      S = mean(S, na.rm = TRUE)
    )

  # Normalizar las frecuencias para cada año
  frecuencias_promedio <- frecuencias_promedio %>%
    mutate(Total = R3 + R2 + R1 + S) %>%
    mutate(across(c(R3, R2, R1, S), ~ ./Total))

  # Transformar e invertir los datos
  datos_largos <- frecuencias_promedio %>%
    pivot_longer(cols = c(S, R1, R2, R3),
                 names_to = "Haplotipo",
                 values_to = "Frecuencia") %>%
    select(-Total) # Eliminar la columna Total, ya que no se necesita más

  return(datos_largos)
}

#devtools::load_all()

#devtools::document()

