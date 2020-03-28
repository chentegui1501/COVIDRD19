a

library(deSolve)

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
      # derivadas parciales
      dS = (-beta * S * I)
      dE = (beta * S * I) - (delta * E)
      dI = (delta * E) - (gamma * I)
      dR = (gamma * I)
      
      # resultados 
      results = c (dS, dE, dI, dR)
      list (results)
    }
  )
}

#Parametros

tasa_contacto = 2.5           # numero de contacto por dias
prob_transmision = .135     # probabilidad de transmision 
period_laten = 14            # periodo latente
period_incu = 7             # periodo incubacion del virus 

beta_value = tasa_contacto * prob_transmision
gamma_value = 1 / period_laten
delta_value = 1 / period_incu

#calcular el r0
R0 = beta_value / gamma_value
R0


#parametros dinamicos de covid
parameter_list = c (beta = beta_value, gamma = gamma_value, delta = delta_value)


beta_value = tasa_contacto * prob_transmision
gamma_value = 1/period_laten
delta_value = 1/period_incu



#calcular el r0


Ro = beta_value / gamma_value



#parametros dinamicos de covid


parametros = c (beta = beta_value, gamma = gamma_value, delta = delta_value)



#valores iniciales 


W = 10000001    # suceptibles_0
X = 1           # infectados_0
Y = 1           # recuperados_0
Z = 0           # expuestos_0

N = W + X + Y + Z #Poblacion



valores_ini = c (S = W/N, E = X/N, I = Y/N, R = Z/N)



tiempo_disc = seq (0, 300, by=1)


output = ode(valores_ini, tiempo_disc, modelo_seir, parametros)



#cambiar out a un data.frame
output <- as.data.frame(output*N) #aqui puede multiplicar 'out' por N
#eliminar la variable 'time' en out
output$time <- NULL



#gráfica
matplot(x = tiempo_disc, y = output, type = "l",
        xlab = "Tiempo", ylab = "S, E, I, R", main = "Modelo SEIR",
        lwd = 1, lty = 1, bty = "l", col = 2:4)
#añadir leyenda 
legend(40, 0.7, c("Susceptibles", "Expuestos", "Infectados", "Recuperados"), 
       pch = 1, col = 2:4, bty = "n", cex = 1)
