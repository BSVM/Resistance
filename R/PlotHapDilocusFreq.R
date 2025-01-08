#' Grafica la frecuencia haplotípica de mutaciones dilocus por año.
#'
#' @param data Datos de la frecuencias haplotípica de mutaciones dilocus a graficar.
#' @param pais País o localidad que se desea graficar.
#'
#' @return Gráfico con la frecuencia haplotípica de mutaciones dilocus por año.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#'
#' @examples
#' PlotHapDilocusFreq(V1016I.F1534C, "Mexico")
#'
#' @export
PlotHapDilocusFreq <- function(data, pais) {
  # Filtrar y transformar los datos
  datos_haplotipos <- Hap_DilocusfreqforYear(data, pais)

  # Reordenar los niveles de Haplotipo
  datos_haplotipos <- datos_haplotipos %>%
    mutate(Haplotipo = factor(Haplotipo, levels = c("R3", "R2", "R1", "S")))

  # Crear el gráfico de barras apiladas con colores personalizados y mejoras en las letras de los ejes
  grafico <- ggplot(datos_haplotipos, aes(x = `Year of Sample`, y = Frecuencia, fill = Haplotipo)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("S" = "#1E90FF", "R1" = "yellow", "R2" = "red", "R3" = "darkblue")) +
    labs(x = "Year", y = "Frequency Haplotipes", fill = "Haplotipe") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 12, face = "plain", angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12, face = "plain", angle = 0, hjust = 1),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold", angle = 90, hjust = 0.5)
    )

  return(grafico)
}

#devtools::load_all()

#devtools::document()

