
datos <- readxl::read_excel("TABLA_ACCIDENTES_20.xlsx")
datos

getwd()

library(ggplot2)
library(dplyr)
library(FSA)

#-------------------------------------------------------------------------------------------------------------------
#########################
### CÓDIGO DE TRABAJO ###
#########################

### Análisis descriptivo de la situación: histogramas iniciales ----

mes <- c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic")
dias_semana <- c("Lun", "Mar", "Mie", "Jue", "Vie", "Sab", "Dom")

# Frecuencia de accidentes en función de la hora del día

ggplot(data = datos, aes(x = HORA)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "cyan", color = "blue", alpha = 0.6) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(breaks = seq(0, 23, by = 1)) +
  labs(title = "nº de accidentes en función de la hora del día", x = "Hora del día", y = "nº de accidentes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

# Frecuencia de accidentes en función del día de la semana

ggplot(data = datos, aes(x = DIA_SEMANA)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "cyan", color = "blue", alpha = 0.6) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(breaks = 1:7, labels = dias_semana)+
  labs(title = "nº de accidentes en función del día de la semana", x = "Día de la semana", y = "nº de accidentes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

# Frecuencia de accidentes en función del mes

ggplot(data = datos, aes(x = MES)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "cyan", color = "blue", alpha = 0.4) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(breaks = 1:12, labels = mes)+
  labs(title = "nº de accidentes en función del mes del año", x = "Mes", y = "nº de accidentes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

# Resumen de parámetros descriptivos

resumen_personalizado <- function(x) {
  res <- c(summary(x), 
           "Desv. estándar" = sd(x, na.rm = TRUE), 
           "IQR" = IQR(x, na.rm = TRUE))
  return(res)
}

resumen_personalizado(datos$HORA); resumen_personalizado(datos$DIA_SEMANA); resumen_personalizado(datos$MES)

# Proporciones de accidentes graves (más de 1 víctima)

lista_1victimaomenos_mes = rep(0, 12)
lista_masde1victima_mes = rep(0, 12)
lista_1victimaomenos_hora = rep(0, 24)
lista_masde1victima_hora = rep(0, 24)
lista_1victimaomenos_diasemana = rep(0, 7)
lista_masde1victima_diasemana = rep(0, 7)

for(i in 1:nrow(datos)) {
  if (datos$TOTAL_VICTIMAS_24H[i] < 2) {lista_1victimaomenos_mes[datos$MES[i]] = lista_1victimaomenos_mes[datos$MES[i]] + 1}
  else {lista_masde1victima_mes[datos$MES[i]] = lista_masde1victima_mes[datos$MES[i]] + 1}
}
for(i in 1:nrow(datos)) {
  if (datos$TOTAL_VICTIMAS_24H[i] < 2) {lista_1victimaomenos_diasemana[datos$DIA_SEMANA[i]] = lista_1victimaomenos_diasemana[datos$DIA_SEMANA[i]] + 1}
  else {lista_masde1victima_diasemana[datos$DIA_SEMANA[i]] = lista_masde1victima_diasemana[datos$DIA_SEMANA[i]] + 1}
}
for(i in 1:nrow(datos)) {
  if (datos$TOTAL_VICTIMAS_24H[i] < 2) {lista_1victimaomenos_hora[datos$HORA[i]+1] = lista_1victimaomenos_hora[datos$HORA[i]+1] + 1}
  else {lista_masde1victima_hora[datos$HORA[i]+1] = lista_masde1victima_hora[datos$HORA[i]+1] + 1}
}
lista_masde1victima_hora
lista_1victimaomenos_hora
proporcion_masde1victima_mes <- lista_masde1victima_mes / lista_1victimaomenos_mes
proporcion_masde1victima_mes
proporcion_masde1victima_hora <- lista_masde1victima_hora / lista_1victimaomenos_hora
proporcion_masde1victima_diasemana <- lista_masde1victima_diasemana / lista_1victimaomenos_diasemana



# Creamos data frame con ambas categorías
df_proporcionmasde1victima_diasemana <- data.frame(
  dias_semana = rep(dias_semana, each = 2),
  categoria = rep(c("Sí", "No"), times = length(dias_semana)),
  proporcion = c(rbind(proporcion_masde1victima_diasemana,
                       1 - proporcion_masde1victima_diasemana))
)
df_proporcionmasde1victima_hora <- data.frame(
  horas = rep(seq(0,23), each = 2),
  categoria = rep(c("Sí", "No"), times = 24),
  proporcion = c(rbind(proporcion_masde1victima_hora, 1 - proporcion_masde1victima_hora))
)
df_proporcionmasde1victima_mes <- data.frame(
  mes = rep(mes, each = 2),
  categoria = rep(c("Sí", "No"), times = length(mes)),
  proporcion = c(rbind(proporcion_masde1victima_mes, 1 - proporcion_masde1victima_mes))
)
df_proporcionmasde1victima_diasemana$dias_semana <- factor(df_proporcionmasde1victima_diasemana$dias_semana, 
                                                           levels = dias_semana)
df_proporcionmasde1victima_hora
df_proporcionmasde1victima_mes$mes <- factor(df_proporcionmasde1victima_mes$mes, levels = mes)

ggplot(df_proporcionmasde1victima_mes, aes(x = mes, y = proporcion, fill = categoria)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#F8766D", "Sí" = "#00BFC4")) +
  labs(
    title = "Accidentes con más de una víctima según mes del año",
    x = "Mes del año",
    y = "Proporción",
    fill = "Más de una víctima"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

ggplot(df_proporcionmasde1victima_diasemana, aes(x = dias_semana, y = proporcion, fill = categoria)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#F8766D", "Sí" = "#00BFC4")) +
  labs(
    title = "Accidentes con más de una víctima según día de la semana",
    x = "Día de la semana",
    y = "Proporción",
    fill = "Más de una víctima"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 17))  

ggplot(df_proporcionmasde1victima_hora, aes(x = horas, y = proporcion, fill = categoria)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("No" = "#F8766D", "Sí" = "#00BFC4")) +
  labs(
    title = "Accidentes con más de una víctima según hora del día",
    x = "Hora del día",
    y = "Proporción",
    fill = "Más de una víctima"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))  

# Test chi cuadrado de Pearson para comprobar la uniformidad de las 
# distribuciones de frecuencia

table_mes <- table(datos$MES)
table_mes
chisq.test(table_mes)

table_hora <- table(datos$HORA)
table_hora
chisq.test(table_hora)

intervalos_diasemana <- cut(datos$DIA_SEMANA, breaks = seq(0.5, 7.5, by = 1), include.lowest = TRUE)
table_diasemana <- as.data.frame(table(Intervalo = intervalos_diasemana))
chisq.test(table_diasemana$Freq)

# Plots de nº total de víctimas en las primeras 24 horas y nº de vehículos involucrados

ggplot(data = datos, aes(x = TOTAL_VICTIMAS_24H)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "cyan", color = "blue", alpha = 0.4) +
  scale_y_continuous(n.breaks = 10) +
  labs(title = "nº de accidentes en función del nº de víctimas en las primeras 24 horas", x = "nº de víctimas", 
       y = "nº de accidentes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

ggplot(data = datos, aes(x = TOTAL_VEHICULOS)) +
  geom_histogram(aes(y = ..count..), binwidth = 1, fill = "cyan", color = "blue", alpha = 0.4) +
  scale_y_continuous(n.breaks = 10) +
  labs(title = "nº de accidentes en función del nº de vehículos involucrados", x = "nº de vehículos involucrados",
       y = "nº de accidentes") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 18))

# Test de Kruskal-Wallis para comprobar si las distribuciones del nº
# de víctimas en las primeras 24 horas y el de vehículos involucrados es 
# igual para los diferentes valores que toman las variables temporales

kruskal.test(TOTAL_VICTIMAS_24H ~ HORA, data = datos)
kruskal.test(TOTAL_VICTIMAS_24H ~ DIA_SEMANA, data = datos)
kruskal.test(TOTAL_VICTIMAS_24H ~ MES, data = datos)

kruskal.test(TOTAL_VEHICULOS ~ HORA, data = datos)
kruskal.test(TOTAL_VEHICULOS ~ DIA_SEMANA, data = datos)
kruskal.test(TOTAL_VEHICULOS ~ MES, data = datos)

# Test post-hoc (Dunn) para extraer más conclusiones

hora_factor <- as.factor(datos$HORA)
dia_semana_factor <- as.factor(datos$DIA_SEMANA)
mes_factor <- as.factor(datos$MES)

dunnTest(TOTAL_VICTIMAS_24H ~ hora_factor, data = datos, method = 'holm')
dunnTest(TOTAL_VEHICULOS ~ hora_factor, data = datos, method = 'holm')
dunnTest(TOTAL_VICTIMAS_24H ~ diasemana_factor, data = datos, method = 'holm')
dunnTest(TOTAL_VEHICULOS ~ diasemana_factor, data = datos, method = 'holm')
dunnTest(TOTAL_VICTIMAS_24H ~ mes_factor, data = datos, method = 'holm')
dunnTest(TOTAL_VEHICULOS ~ mes_factor, data = datos, method = 'holm')

dunntest_diasemana_24h <- dunnTest(TOTAL_VICTIMAS_24H ~ diasemana_factor, 
                                   data = datos, method = 'holm')$res
dunntest_diasemana_vehiculos <- dunnTest(TOTAL_VEHICULOS ~ diasemana_factor, 
                                         data = datos, method = 'holm')$res

dunntest_diasemana_24h$Comparison <- factor(dunntest_diasemana_24h$Comparison, 
                                            levels = unique(dunntest_diasemana_24h$Comparison))
dunntest_diasemana_24h
dunntest_diasemana_vehiculos$Comparison <- factor(dunntest_diasemana_vehiculos$Comparison, 
                                                  levels = unique(dunntest_diasemana_vehiculos$Comparison))
dunntest_diasemana_vehiculos

# Se aprecia un p-valor mucho mayor al comparar las distribuciones de días laborales. De 
# forma más visual, se representan los p-valores ajustados frente al nivel de significación

ggplot(dunntest_diasemana_24h, aes(x = Comparison, y = P.adj)) +
  geom_point(color = "blue", size = 4) +
  scale_y_log10() +
  labs(title = "comparación por pares del nº de víctimas en 24 horas según el día de la semana",
       x = "Pares comparados",
       y = "p-valores (escala logarítmica)") +
  geom_hline(aes(yintercept = 0.01),
             color = "black",
             linetype = "dashed", size = 1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))

ggplot(dunntest_diasemana_vehiculos, aes(x = Comparison, y = P.adj)) +
  geom_point(color = "blue", size = 4) +
  scale_y_log10() +
  labs(title = "comparación por pares del nº de vehículos involucrados según el día de la semana",
       x = "Pares comparados",
       y = "p-valores (escala logarítmica)") +
  geom_hline(aes(yintercept = 0.01),
             color = "black",
             linetype = "dashed", size = 1) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16))

# Se vuelve a aplicar Kruskal-Wallis, pero analizando los días laborales y el 
# par sábado domingo por separado

datos_dias_laborales <- datos[datos$DIA_SEMANA %in% c(1,2,3,4,5),]
datos_findesemana <- datos[datos$DIA_SEMANA %in% c(6,7),]

kruskal.test(TOTAL_VEHICULOS ~ DIA_SEMANA, data = datos_dias_laborales)
kruskal.test(TOTAL_VICTIMAS_24H ~ DIA_SEMANA, data = datos_dias_laborales)

kruskal.test(TOTAL_VEHICULOS ~ DIA_SEMANA, data = datos_findesemana)
kruskal.test(TOTAL_VICTIMAS_24H ~ DIA_SEMANA, data = datos_findesemana)