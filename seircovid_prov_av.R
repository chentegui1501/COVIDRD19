
########## Codigo: Modelo SEIR COVID19: Republica Dominicana, con Analisis de Politicas de Movibildad######
#########  Francisco A. Ramírez / Vicente Peña Peralta
#########  Universidad Autónoma de Santo Domingo
rm(list=ls())
#Librerias a utilizar
library(deSolve)
library(ggplot2)
library(tseries)
library(forecast)
library(readxl)

d    <- ts(read_excel("data_para_politicas.xlsx"),freq =365, start = c(2020,62))


#Modelo SEIR

modelo_seir = function (valores_iniciales, valores_estado, parametros)
{
  # creamos las variables de estados
  S_DN = valores_estado [1]        # susceptibles
  S_DU = valores_estado [2]
  S_SD = valores_estado [3]
  S_ST = valores_estado [4]
  S_RE = valores_estado [5]
  
  E_DN = valores_estado [6]
  E_DU = valores_estado [7]
  E_SD = valores_estado [8]
  E_ST = valores_estado [9]
  E_RE = valores_estado [10]
  
  I_DN = valores_estado [11]
  I_DU = valores_estado [12]
  I_SD = valores_estado [13]
  I_ST = valores_estado [14]
  I_RE = valores_estado [15]
  
  R_DN = valores_estado [16]        # recuperados
  R_DU = valores_estado [17]        # recuperados
  R_SD = valores_estado [18]        # recuperados
  R_ST = valores_estado [19]        # recuperados
  R_RE = valores_estado [20]        # recuperados
  
  S    = valores_estado [21]
  E    = valores_estado [22]
  I    = valores_estado [23]
  R    = valores_estado [24]
  
  with ( 
    as.list (parametros), 
    {
      b_dn = b0_dn*(1-a0_dn)*(1-0.05*I_DN)^k0_dn
      b_sd = b0_sd*(1-a0_sd)*(1-0.05*I_SD)^k0_sd
      b_du = b0_du*(1-a0_du)*(1-0.05*I_DU)^k0_du
      b_st = b0_st*(1-a0_st)*(1-0.05*I_ST)^k0_st
      b_re = b0_re*(1-a0_re)*(1-0.05*I_RE)^k0_re
      
      # derivadas parciales
      dS_DN = (-b_dn * S_DN * I_DN)
      dS_DU = (-b_du * S_DU * I_DU)
      dS_SD = (-b_sd * S_SD * I_SD)
      dS_ST = (-b_st * S_ST * I_ST)
      dS_RE = (-b_re * S_RE * I_RE)
      
      dE_DN = (b_dn * S_DN * I_DN) - (delta * E_DN)
      dE_DU = (b_du * S_DU * I_DU) - (delta * E_DU)
      dE_SD = (b_sd * S_SD * I_SD) - (delta * E_SD)
      dE_ST = (b_st * S_ST * I_ST) - (delta * E_ST)
      dE_RE = (b_re * S_RE * I_RE) - (delta * E_RE)
      
      dI_DN = (delta * E_DN) - (gamma * I_DN)
      dI_DU = (delta * E_DU) - (gamma * I_DU)
      dI_SD = (delta * E_SD) - (gamma * I_SD)
      dI_ST = (delta * E_ST) - (gamma * I_ST)
      dI_RE = (delta * E_RE) - (gamma * I_RE)
      
      dR_DN = (gamma * I_DN)
      dR_DU = (gamma * I_DU)
      dR_SD = (gamma * I_SD)
      dR_ST = (gamma * I_ST)
      dR_RE = (gamma * I_RE)
      
      dS = dS_DN + dS_DU + dS_SD + dS_ST + dS_RE
      dE = dE_DN + dE_DU + dE_SD + dE_ST + dE_RE
      dI = dI_DN + dI_DU + dI_SD + dI_ST + dI_RE
      dR = dR_DN + dR_DU + dR_SD + dR_ST + dR_RE
      
      # resultados 
      results = c (dS_DN, dS_DU, dS_SD, dS_ST, dS_RE,
                   dE_DN, dE_DU, dE_SD, dE_ST, dE_RE, 
                   dI_DN, dI_DU, dI_SD, dI_ST, dI_RE, 
                   dR_DN, dR_DU, dR_SD, dR_ST, dR_RE,
                   dS, dE, dI, dR)
      list (results)
    }
  )
}



#Parametros Calibrados Generales
period_laten <- 14            # periodo latente
period_incu  <- 7              # periodo incubacion del virus 
gamma_value  <- 1 / period_laten
delta_value  <- 1 / period_incu
tiempo_disc <- seq (0, 300, by=1)

a0_dn  <- 0
a0_du  <- 0
a0_sd  <- 0
a0_st  <- 0
a0_re  <- 0

b0_dn  <- 0.85 
b0_du  <- 0.58
b0_sd  <- 0.61
b0_st  <- 0.62
b0_re  <- 0.72
  
k0_dn  <- 100
k0_du  <- 100
k0_sd  <- 100
k0_st  <- 100
k0_re  <- 100

p <- c(a0_dn = a0_dn, a0_du = a0_du, a0_sd = a0_sd, a0_st = a0_st, a0_re = a0_re,
       b0_dn = b0_dn, b0_du = b0_du, b0_sd = b0_sd, b0_st = b0_st, b0_re = b0_re,
       k0_dn = k0_dn, k0_du = k0_du, k0_sd = k0_sd, k0_st = k0_st, k0_re = k0_re,
       gamma = gamma_value, delta = delta_value)
#Suceptibles
W_DN <- 965040
W_DU <- 289574
W_SD <- 2374370
W_ST <- 963422
W_RE <- 4852875
W    <- W_DN + W_DU + W_SD + W_ST + W_RE
#Expuestps
X_DN <- 0
X_DU <- 0
X_SD <- 0
X_ST <- 0
X_RE <- 0
X    <- X_DN + X_DU + X_SD + X_ST + X_RE 
#Infectados
Y_DN <- 1
Y_DU <- 1
Y_SD <- 1
Y_ST <- 1
Y_RE <- 1
Y    <- Y_DN + Y_DU + Y_SD + Y_ST + Y_RE
#Recuperados
Z_DN <- 0
Z_DU <- 0
Z_SD <- 0
Z_ST <- 0
Z_RE <- 0
Z   <- Z_DN + Z_DU + Z_SD + Z_ST +Z_RE
#Poblaciones
N_DN <- W_DN + X_DN + Y_DN + Z_DN
N_DU <- W_DU + X_DU + Y_DU + Z_DU
N_SD <- W_SD + X_SD + Y_SD + Z_SD
N_ST <- W_ST + X_ST + Y_ST + Z_ST
N_RE <- W_RE + X_RE + Y_RE + Z_RE 
N <- N_DN + N_DU + N_SD + N_ST + N_RE

valor_ini <- c(S_DN = W_DN/N_DN, S_DU = W_DU/N_DU, S_SD = W_SD/N_SD, S_ST = W_ST/N_ST, S_RE = W_RE/N_RE,
               E_DN = X_DN/N_DN, E_DU = X_DU/N_DU, E_SD = X_SD/N_SD, E_ST = X_ST/N_ST, E_RE = X_RE/N_RE,
               I_DN = Y_DN/N_DN, I_DU = Y_DU/N_DU, I_SD = Y_SD/N_SD, I_ST = Y_ST/N_ST, I_RE = Y_RE/N_RE,
               R_DN = Z_DN/N_DN, R_DU = Z_DU/N_DU, R_SD = Z_SD/N_SD, R_ST = Z_ST/N_ST, R_RE = Z_RE/N_RE, 
               S = W/N, E = X/N, I = Y/N, R = Z/N)

#valor_ini <- c(S_DN = W_DN/N, S_DU = W_DU/N, S_SD = W_SD/N, S_ST = W_ST/N, S_RE = W_RE/N,
#               E_DN = X_DN/N, E_DU = X_DU/N, E_SD = X_SD/N, E_ST = X_ST/N, E_RE = X_RE/N,
#               I_DN = Y_DN/N, I_DU = Y_DU/N, I_SD = Y_SD/N, I_ST = Y_ST/N, I_RE = Y_RE/N,
#               R_DN = Z_DN/N, R_DU = Z_DU/N, R_SD = Z_SD/N, R_ST = Z_ST/N, R_RE = Z_RE/N)


#Simulacion
o  <-  ode(valor_ini, tiempo_disc, modelo_seir, p)
o  <- as.data.frame(o)
o$time <- NULL
st  <- ts(o, freq = 365, start = c(2020,62))
         


#----------------
#Distrito Nacional
#----------------

matplot(x = tiempo_disc, y = st[,c(1,6,11,16)]*N_DN, type = "l",
        xlab = "Días", ylab = "SEIR", main = "Modelo SEIR, Distrito Nacional",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(100, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)

#Ajuste de los Casos Publicados por Salud Publica
s_dn <- ts(st[1:30,11]*N_DN,freq = 365, start = c(2020,62))
d_dn <- d[,2]

data_dn <-  window(cbind(d_dn,s_dn),start=c(2020,80))

data_dn %>% autoplot(xlab = "Dias",
                     ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo, Distrito Nacional")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()



#----------------
#Duarte
#----------------

matplot(x = tiempo_disc, y = st[,c(2,7,12,17)]*N_DU, type = "l",
        xlab = "Días", ylab = "SEIR", main = "Modelo SEIR, Duarte",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(100, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)

#Ajuste de los Casos Publicados por Salud Publica
s_du <- ts(st[1:30,12]*N_DU,freq = 365, start = c(2020,62))
d_du <- d[,3]

data_du <-  window(cbind(d_du,s_du),start=c(2020,80))

data_du %>% autoplot(xlab = "Dias",
                     ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo, Duarte")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()


#----------------
#Santo Domingo
#----------------

matplot(x = tiempo_disc, y = st[,c(3,8,13,18)]*N_SD, type = "l",
        xlab = "Días", ylab = "SEIR", main = "Modelo SEIR, Santo Domingo",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(100, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)

#Ajuste de los Casos Publicados por Salud Publica
s_sd <- ts(st[1:30,13]*N_SD,freq = 365, start = c(2020,62))
d_sd <- d[,4]

data_sd <-  window(cbind(d_sd,s_sd),start=c(2020,80))

data_sd %>% autoplot(xlab = "Dias",
                     ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo, Santo Domingo")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()


#----------------
#Santiago
#----------------

matplot(x = tiempo_disc, y = st[,c(4,9,14,19)]*N_ST, type = "l",
        xlab = "Días", ylab = "SEIR", main = "Modelo SEIR, Santiago",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(100, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)

#Ajuste de los Casos Publicados por Salud Publica
s_st <- ts(st[1:30,14]*N_ST,freq = 365, start = c(2020,62))
d_st <- d[,5]

data_st <-  window(cbind(d_st,s_st),start=c(2020,80))

data_st %>% autoplot(xlab = "Dias",
                     ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo, Santiago")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()


#Resto del Pais
#----------------

matplot(x = tiempo_disc, y = st[,c(5,10,15,20)]*N_RE, type = "l",
        xlab = "Días", ylab = "SEIR", main = "Modelo SEIR, Resto del País",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(100, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)

#Ajuste de los Casos Publicados por Salud Publica
s_re <- ts(st[1:30,15]*N_RE,freq = 365, start = c(2020,62))
d_re <- d[,6]

data_re <-  window(cbind(d_re,s_re),start=c(2020,80))

data_re %>% autoplot(xlab = "Dias",
                     ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo, Resto del País")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()

#Graficando
#--------------
#Pais
#-------------

#Ajuste de los Casos Publicados por Salud Publica
s <- s_dn+s_du+s_re+s_sd+s_st
d_ <- d[,1]

data_ <-  window(cbind(d_,s),start=c(2020,80))

data_ %>% autoplot(xlab = "Dias",
                   ylab = "Infectados")+
  ggtitle(label = "Evolución de los Casos de Infectados",
          subtitle = "Desempeño del Modelo")+
  labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
  theme_bw()+
  geom_point()



## Analisis de las politicas
#Suceptibles
W_DN <- 965040
W_DU <- 289574
W_SD <- 2374370
W_ST <- 963422
W_RE <- 4852875
W    <- W_DN + W_DU + W_SD + W_ST + W_RE
#Expuestps
X_DN <- 0
X_DU <- 0
X_SD <- 0
X_ST <- 0
X_RE <- 0
X    <- X_DN + X_DU + X_SD + X_ST + X_RE 
#Infectados
Y_DN <- 462
Y_DU <- 94
Y_SD <- 145
Y_ST <- 122
Y_RE <- 286
Y    <- Y_DN + Y_DU + Y_SD + Y_ST + Y_RE
#Recuperados
Z_DN <- 0
Z_DU <- 0
Z_SD <- 0
Z_ST <- 0
Z_RE <- 0
Z   <- Z_DN + Z_DU + Z_SD + Z_ST +Z_RE
#Poblaciones
N_DN <- W_DN + X_DN + Y_DN + Z_DN
N_DU <- W_DU + X_DU + Y_DU + Z_DU
N_SD <- W_SD + X_SD + Y_SD + Z_SD
N_ST <- W_ST + X_ST + Y_ST + Z_ST
N_RE <- W_RE + X_RE + Y_RE + Z_RE 
N <- N_DN + N_DU + N_SD + N_ST + N_RE

valor_ini <- c(S_DN = W_DN/N_DN, S_DU = W_DU/N_DU, S_SD = W_SD/N_SD, S_ST = W_ST/N_ST, S_RE = W_RE/N_RE,
               E_DN = X_DN/N_DN, E_DU = X_DU/N_DU, E_SD = X_SD/N_SD, E_ST = X_ST/N_ST, E_RE = X_RE/N_RE,
               I_DN = Y_DN/N_DN, I_DU = Y_DU/N_DU, I_SD = Y_SD/N_SD, I_ST = Y_ST/N_ST, I_RE = Y_RE/N_RE,
               R_DN = Z_DN/N_DN, R_DU = Z_DU/N_DU, R_SD = Z_SD/N_SD, R_ST = Z_ST/N_ST, R_RE = Z_RE/N_RE, 
               S = W/N, E = X/N, I = Y/N, R = Z/N)


### Escenario 1: Cuarentena total en DN y Duarte, y parcial resto del pais
a0_dn  <- 0.8
a0_du  <- 0.8
a0_sd  <- 0.4
a0_st  <- 0.4
a0_re  <- 0.4

p_e1 <- c(a0_dn = a0_dn, a0_du = a0_du, a0_sd = a0_sd, a0_st = a0_st, a0_re = a0_re,
       b0_dn = b0_dn, b0_du = b0_du, b0_sd = b0_sd, b0_st = b0_st, b0_re = b0_re,
       k0_dn = k0_dn, k0_du = k0_du, k0_sd = k0_sd, k0_st = k0_st, k0_re = k0_re,
       gamma = gamma_value, delta = delta_value)

oe1  <-  ode(valor_ini, tiempo_disc, modelo_seir, p_e1)
oe1  <- as.data.frame(oe1)
oe1$time <- NULL

st_e1  <- ts(oe1, freq = 365, start = c(2020,62))

s_dn_e1 <- ts(st_e1[1:15,11]*N_DN,freq = 365, start = c(2020,62))
s_du_e1 <- ts(st_e1[1:15,12]*N_DU,freq = 365, start = c(2020,62))
s_sd_e1 <- ts(st_e1[1:15,13]*N_SD,freq = 365, start = c(2020,62))
s_st_e1 <- ts(st_e1[1:15,14]*N_ST,freq = 365, start = c(2020,62))
s_re_e1 <- ts(st_e1[1:15,15]*N_RE,freq = 365, start = c(2020,62))
s_e1 <- s_dn_e1 + s_du_e1 + s_re_e1 + s_sd_e1 + s_st_e1



### Escenario 2: Cuarentena todo el pais

a0_dn  <- 0.8
a0_du  <- 0.8
a0_sd  <- 0.8
a0_st  <- 0.8
a0_re  <- 0.8

p_e2 <- c(a0_dn = a0_dn, a0_du = a0_du, a0_sd = a0_sd, a0_st = a0_st, a0_re = a0_re,
          b0_dn = b0_dn, b0_du = b0_du, b0_sd = b0_sd, b0_st = b0_st, b0_re = b0_re,
          k0_dn = k0_dn, k0_du = k0_du, k0_sd = k0_sd, k0_st = k0_st, k0_re = k0_re,
          gamma = gamma_value, delta = delta_value)

oe2  <-  ode(valor_ini, tiempo_disc, modelo_seir, p_e2)
oe2  <- as.data.frame(oe2)
oe2$time <- NULL
st_e2  <- ts(oe2, freq = 365, start = c(2020,62))

s_dn_e2 <- ts(st_e2[1:15,11]*N_DN,freq = 365, start = c(2020,62))
s_du_e2 <- ts(st_e2[1:15,12]*N_DU,freq = 365, start = c(2020,62))
s_sd_e2 <- ts(st_e2[1:15,13]*N_SD,freq = 365, start = c(2020,62))
s_st_e2 <- ts(st_e2[1:15,14]*N_ST,freq = 365, start = c(2020,62))
s_re_e2 <- ts(st_e2[1:15,15]*N_RE,freq = 365, start = c(2020,62))
s_e2 <- s_dn_e2 + s_du_e2 + s_re_e2 + s_sd_e2 + s_st_e2



### Graficando

#data  <- ts(read_excel("datos_salud_publica.xlsx"),freq =365, start = c(2020,62))
startna <- ts(rep(NA,13), start = c(2020,62))
esc1 <- ts(c(d[18:30,1],s_e1[1:15]),freq = 365, start = c(2020,62))
esc2 <- ts(c(startna,s_e2[1:15]),freq = 365, start = c(2020,62))


autoplot(cbind(esc1,esc2),
         xlab = "Días desde la medida",
         ylab = "Infectados")+
  ggtitle(label = "Impacto de las Medidas de Movibilidad sobre Dinámica de los Infectados",
          subtitle = "Simulaciones Modelo SEIR")+
  geom_point(size =1)+
  #labs(caption = "Datos de los Boletines de Salud Pública")+
  scale_color_manual( name = "Medidas de Distanciamiento", label = c("Escenario 1", "Escenario 2"), 
                      values = c("black","red","green", "blue"))+
  theme_bw()+
  theme(legend.position = "top")

