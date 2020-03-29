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

d    <- ts(read_excel("datos_sp_provincias.xlsx"),freq =365, start = c(2020,62))


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


#Parametros Calibrados Generales
period_laten = 14            # periodo latente
period_incu = 7             # periodo incubacion del virus 
gamma_value = 1 / period_laten
delta_value = 1 / period_incu
R0 = beta_value / gamma_value   #R0 Implicito calibración inicial
alpha0 = 0
kappa0 = 100
tiempo_disc = seq (0, 300, by=1)

###----ANALISIS POR PROVINCIAS -----

# Distrito Nacional
  beta_dn <-  .85
  p_dn = c (alpha0 = alpha0, beta0 = beta_dn, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
  #p_dn = c (beta = beta_dn, gamma = gamma_value)
  W_dn = 965040 -1  # suceptibles_0
  X_dn = 0           # Expuestos_0
  Y_dn = 1        # Infectados_0
  Z_dn = 0           # Recuperados_0
  N_dn = W_dn + X_dn + Y_dn + Z_dn #Poblacion

  valores_ini = c (S = W_dn/N_dn, E = X_dn/N_dn, I = Y_dn/N_dn, R = Z_dn/N_dn)
  
  #Solucionando el Sistema 
  o_dn = ode(valores_ini, tiempo_disc, modelo_seir, p_dn)
  o_dn <- as.data.frame(o_dn*N_dn) 

  #Ajuste de los Casos Publicados por Salud Publica
  s_dn <- ts(o_dn[1:28,4],freq = 365, start = c(2020,62))
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Santo Domingo
  beta_sd <-  0.61
  p_sd = c (alpha0 = alpha0, beta0 = beta_sd, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
  W_sd = 2374370 -1  # suceptibles_0
  X_sd = 0           # Expuestos_0
  Y_sd = 1        # Infectados_0
  Z_sd = 0           # Recuperados_0
  N_sd = W_sd + X_sd + Y_sd + Z_sd #Poblacion
  
  valores_ini = c (S = W_sd/N_sd, E = X_sd/N_sd, I = Y_sd/N_sd, R = Z_sd/N_sd)
  
  #Solucionando el Sistema 
  o_sd = ode(valores_ini, tiempo_disc, modelo_seir, p_sd)
  o_sd <- as.data.frame(o_sd*N_sd) 
  
  #Ajuste de los Casos Publicados por Salud Publica
  s_sd <- ts(o_sd[1:28,4],freq = 365, start = c(2020,62))
  d_sd <- d[,33]
  
  data_sd <- window(cbind(d_sd,s_sd),start=c(2020,80))
  
  data_sd %>% autoplot(xlab = "Dias",
                       ylab = "Infectados")+
    ggtitle(label = "Evolución de los Casos de Infectados",
            subtitle = "Desempeño del Modelo, Santo Domingo")+
    labs(caption = "Datos de los Boletines de Salud Pública")+
    scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
    theme_bw()+
    geom_point()
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Duarte
  beta_du <-  0.58
  p_du = c (alpha0 = alpha0, beta0 = beta_du, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
  W_du = 289574 -1  # suceptibles_0
  X_du = 0           # Expuestos_0
  Y_du = 1         # Infectados_0
  Z_du = 0           # Recuperados_0
  N_du = W_du + X_du + Y_du + Z_du #Poblacion
  
  valores_ini = c (S = W_du/N_du, E = X_du/N_du, I = Y_du/N_du, R = Z_du/N_du)
  
  #Solucionando el Sistema 
  o_du = ode(valores_ini, tiempo_disc, modelo_seir, p_du)
  o_du <- as.data.frame(o_du*N_du) 
  
  #Ajuste de los Casos Publicados por Salud Publica
  s_du <- ts(o_du[1:28,4],freq = 365, start = c(2020,62))
  d_du <- d[,7]
  
  data_du <- window(cbind(d_du,s_du),start=c(2020,80))
  
  data_du %>% autoplot(xlab = "Dias",
                       ylab = "Infectados")+
    ggtitle(label = "Evolución de los Casos de Infectados",
            subtitle = "Desempeño del Modelo, Duarte")+
    labs(caption = "Datos de los Boletines de Salud Pública")+
    scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
    theme_bw()+
    geom_point()
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Santiago
  beta_st <-  0.62
  p_st = c (alpha0 = alpha0, beta0 = beta_st, kappa0 = kappa0, gamma = gamma_value, delta = delta_value)
  W_st = 963422 -1  # suceptibles_0
  X_st = 0           # Expuestos_0
  Y_st = 1         # Infectados_0
  Z_st = 0           # Recuperados_0
  N_st = W_st + X_st + Y_st + Z_st #Poblacion
  
  valores_ini = c (S = W_st/N_st, E = X_st/N_st, I = Y_st/N_st, R = Z_st/N_st)
  
  #Solucionando el Sistema 
  o_st = ode(valores_ini, tiempo_disc, modelo_seir, p_st)
  o_st <- as.data.frame(o_st*N_st) 
  
  #Ajuste de los Casos Publicados por Salud Publica
  s_st <- ts(o_st[1:28,4],freq = 365, start = c(2020,62))
  d_st <- d[,26]
  
  data_st <- window(cbind(d_st,s_st),start=c(2020,80))
  
  data_st %>% autoplot(xlab = "Dias",
                       ylab = "Infectados")+
    ggtitle(label = "Evolución de los Casos de Infectados",
            subtitle = "Desempeño del Modelo, Santiago")+
    labs(caption = "Datos de los Boletines de Salud Pública")+
    scale_color_manual( name = element_blank(), label = c("Datos", "Modelo"), values = c("blue", "red"))+
    theme_bw()+
    geom_point()
  
  
  
  