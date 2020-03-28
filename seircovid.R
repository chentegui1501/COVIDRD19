########## Codigo: Modelo SEIR COVID19: Republica Dominicana, con Analisis de Politicas de Movibildad######
#########  Francisco A. Ramírez / Vicente Peña Peralta
#########  Universidad Autónoma de Santo Domingo

#Librerias a utilizar
library(deSolve)
library(ggplot2)
library(tseries)
library(forecast)
library(readxl)

#Modelo SEIR

modelo_seir = function (valores_iniciales, valores_estado, parametros)
{
  # creamos las variables de estados
  S = valores_estado [1]        # susceptibles
  E = valores_estado [2]        # expuestos
  I = valores_estado [3]        # infectados
  R = valores_estado [4]        # recuperados
  
  with ( 
    as.list (parametros), 
    {
      beta = beta0*(1-alpha0)*(1-0.05*I)^kappa0
      # derivadas parciales
      dS = (-beta * S * I)
      dE = (beta * S * I) - (delta * E)
      dI = (delta * E) - (gamma * I)
      dR = (gamma * I)
      
      tasa_contag = beta/gamma
      # resultados 
      results = c (dS, dE, dI, dR)
      list (results,tasa_contag)
    }
  )
}

#Parametros Calibrados
period_laten = 14            # periodo latente
period_incu = 7             # periodo incubacion del virus 
beta_value =.92    #0.905
gamma_value = 1 / period_laten
delta_value = 1 / period_incu
R0 = beta_value / gamma_value   #R0 Implicito calibración inicial
alpha0 = 0
kappa0 = 100


parametros = c (alpha0 = alpha0, beta0 = beta_value, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)

#Valores iniciales 


W = 10000000 -1 # suceptibles_0
X = 1           # Expuestos_0
Y = 1           # Infectados_0
Z = 0           # Recuperados_0

N = W + X + Y + Z #Poblacion

valores_ini = c (S = W/N, E = X/N, I = Y/N, R = Z/N)

tiempo_disc = seq (0, 300, by=1)

#Solucionando el Sistema 
output = ode(valores_ini, tiempo_disc, modelo_seir, parametros)

#cambiar out a un data.frame
output <- as.data.frame(output*N) #aqui puede multiplicar 'out' por N
#eliminar la variable 'time' en out
output$time <- NULL



#Gráfica de la Simulación
matplot(x = tiempo_disc, y = output, type = "l",
        xlab = "Tiempo", ylab = "S, E, I, R", main = "Modelo SEIR",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
#añadir leyenda 
legend(40, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)


#Ajuste de los Casos Publicados por Salud Publica

simul <- ts(output[1:26,3],frequency = 365, start = c(2020,62))
data  <- ts(read_excel("datos_salud_publica.xlsx"),freq =365, start = c(2020,62))

autoplot(cbind(data,simul),
         xlab = "Dias",
         ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",subtitle = "Desempeño del Modelo")+
  geom_point()+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()



#Simulaciones de Politica

## Valores Iniciales: Datos al dia 26

W = 9999200  # suceptibles_0
X = 422           # Expuestos_0
Y = 638           # Infectados_0
Z = 158           # Recuperados_0

N = W + X + Y + Z #Poblacion
N = 10000000

valores_ini = c (S = W/N, E = X/N, I = Y/N, R = Z/N)


alpha0 = 0
parametros = c (alpha0 = alpha0, beta0 = beta_value, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
output1 = ode(valores_ini, tiempo_disc, modelo_seir, parametros)

alpha0 = 0.4
parametros = c (alpha0 = alpha0, beta0 = beta_value, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
output2 = ode(valores_ini, tiempo_disc, modelo_seir, parametros)


alpha0 = 0.6
parametros = c (alpha0 = alpha0, beta0 = beta_value, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
output3 = ode(valores_ini, tiempo_disc, modelo_seir, parametros)


alpha0 = 0.8
parametros = c (alpha0 = alpha0, beta0 = beta_value, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
output4 = ode(valores_ini, tiempo_disc, modelo_seir, parametros)


#Grafico 15 dias
#data <- ts(data, freq = 365, start = c(2020, 62), end = c(2020, 102))
startna <- ts(rep(NA,26), start = c(2020,62))
esc1 <- ts(c(data,output1[1:15,4]*N),freq = 365, start = c(2020,62))
esc2 <- ts(c(startna,output2[1:15,4]*N),freq = 365, start = c(2020,62))
esc3 <- ts(c(startna,output3[1:15,4]*N),freq = 365, start = c(2020,62))
esc4 <- ts(c(startna,output4[1:15,4]*N),freq = 365, start = c(2020,62))

simulc <- ts(cbind(esc1,esc2,esc3,esc4),
             freq = 365, start=c(2020,62 ))


autoplot(cbind(simulc),
         xlab = "Días desde la medida",
         ylab = "Infectados")+
  ggtitle(label = "Impacto de las Medidas de Movibilidad sobre Dinámica de los Infectados",
          subtitle = "Simulaciones Modelo SEIR")+
  geom_point()+
  #labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = "Medidas de Distanciamiento", label = c("Ninguna", "Suave", "Duras", "Extremas"), 
                      values = c("black","red","green", "blue"))+
  theme_bw()+
  theme(legend.position = "top")

#Grafico 90 dias

startna <- ts(rep(NA,26), start = c(2020,62))
esc1 <- ts(c(data,output1[1:90,4]*N),freq = 365, start = c(2020,62))
esc2 <- ts(c(startna,output2[1:90,4]*N),freq = 365, start = c(2020,62))
esc3 <- ts(c(startna,output3[1:90,4]*N),freq = 365, start = c(2020,62))
esc4 <- ts(c(startna,output4[1:90,4]*N),freq = 365, start = c(2020,62))

simuld <- ts(cbind(esc1,esc2,esc3,esc4),
             freq = 365, start=c(2020,62 ))


autoplot(cbind(simuld*N),
         xlab = "Días desde la medida",
         ylab = "Infectados")+
  ggtitle(label = "Impacto de las Medidas de Movibilidad sobre Dinámica de los Infectados",
          subtitle = "Simulaciones Modelo SEIR")+
  geom_point()+
  #labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = "Medidas de Distanciamiento", label = c("Ninguna", "Suave", "Duras", "Extremas"), 
                      values = c("black","red","green", "blue"))+
  theme_bw()+
  theme(legend.position = "top")



#Grafico 15 dias de tasa de contagio

startna <- ts(rep(NA,26), start = c(2020,62))
esc1 <- ts(output1[1:15,5],freq = 365, start = c(2020,62))
esc2 <- ts(output2[1:15,5],freq = 365, start = c(2020,62))
esc3 <- ts(output3[1:15,5],freq = 365, start = c(2020,62))
esc4 <- ts(output4[1:15,5],freq = 365, start = c(2020,62))

simulc <- ts(cbind(esc1,esc2,esc3,esc4),
             freq = 365, start=c(2020,62 ))


autoplot(cbind(simulc),
         xlab = "Días desde la medida",
         ylab = "Infectados")+
  ggtitle(label = "Impacto de las Medidas de Movibilidad sobre Dinámica de los Infectados",
          subtitle = "Simulaciones Modelo SEIR")+
  geom_point()+
  #labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = "Medidas de Distanciamiento", label = c("Ninguna", "Suave", "Duras", "Extremas"), 
                      values = c("black","red","green", "blue"))+
  theme_bw()+
  theme(legend.position = "top")


