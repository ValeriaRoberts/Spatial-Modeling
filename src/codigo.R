##################### Carga de librerias #######################

# Proposito general
library(tidyverse)
library(ggplot2)

# Simuladores MCMC
library(R2OpenBUGS)
library(R2jags)

# Graficas espaciales
library(SpatialEpi)
library(maps)
library(maptools)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sp)

#################### Funciones auxiliares ###################### 

prob <- function(x){
  out<-min(length(x[x>0])/length(x),length(x[x<0])/length(x))
  out
}

get_x <- function(value){(value |> strsplit(" "))[[1]][1]}
get_y <- function(value){(value |> strsplit(" "))[[1]][2]}

# La siguiente funcion nos ayuda a procesar las coordenadas
# que vienen en el csv, tal que la libreria sp los procese
# como poligonos
obtener_poligono <- function(pais){
  
  # Filtramos el pais y sacamos el poligono como viene
  poligonos <- datos |>
    filter(Country == pais) |>
    head(1) |>
    pull(geometry) |>
    strsplit("\\(")
  
  ### Como hay algunos multipoligonos, queremos quedarnos con el
  ### mas grande
  
  # Primero vemos las longitudes de cada poligono
  poligonos_length <- c()
  for (poligono in poligonos[[1]]){
    poligonos_length <- c(poligonos_length, poligono |> nchar())
  }
  
  # Al mayor le aplicamos un procesamiento para tener las
  # coordenadas como las necesitamos
  country_polygon <- poligonos[[1]][which.max(poligonos_length)]
  country_polygon <- gsub('\\)', '', country_polygon)
  country_polygon <- country_polygon |> strsplit("\\,")
  country_polygon <- country_polygon[[1]] |> as_tibble()
  country_polygon <- country_polygon |>
    mutate(value = ifelse(value |> substr(1, 1) == " ",
                          value |> substr(2, 1000),
                          value))
  
  # Obtenemos las latitudes (x) y longitudes (y)
  x <- country_polygon  |> pull(value) |> sapply(get_x)
  y <- country_polygon  |> pull(value) |> sapply(get_y)
  
  # Y las guardamos como numericos en una matriz
  country_polygon <- country_polygon |> 
    mutate(x_coord = as.numeric(x), y_coord = as.numeric(y)) |>
    select(-c("value")) |>
    as.matrix()
  
  # Creamos el objeto Polygon
  p = Polygon(country_polygon |> as.matrix())
  ps = Polygons(list(p), pais)
  
  ps
  
}

plot_arrows <- function(region){
  
  if (region == region){
    arrows(-75, -39, -72, -38, length = 0.1) #Chile
    arrows(-56, 8, -59, 5, length = 0.1) #Guyana
    arrows(-81, 2, -78, -1, length = 0.1) #Ecuador
    arrows(-75, -22.5, -59, -22.5, length = 0.1) #Paraguay
    arrows(-52, -35, -56, -32, length = 0.1) #Uruguay
  }
  
}

save_definition <- function(model_id, data, inits, parameters, covid.sim){
  
  dir.create(paste0("../results/",model_id))
  
  sink(paste0("../results/", model_id, "/data.txt"))
  print(data)
  sink()
  
  sink(paste0("../results/", model_id, "/inits.txt"))
  print(inits)
  sink()
  
  sink(paste0("../results/", model_id, "/parameters.txt"))
  print(parameters)
  sink()
  
  sink(paste0("../results/", model_id, "/dic.txt"))
  print(covid.sim$DIC)
  sink()
  
  save_mcmc_params(model_id, covid.sim)
  
  save_traceplots(model_id, covid.sim)
}

save_traceplots <- function(model_id, covid.sim){
  
  # chain <- c()
  chain <- rep(1:covid.sim$n.chains, covid.sim$n.sims / covid.sim$n.chains)
  n_sim <- c()
  
  for (ch in 1:covid.sim$n.chains){
    # chain <- c(chain, rep(ch, covid.sim$n.sims / covid.sim$n.chains))
    n_sim <- c(n_sim, seq(1:(covid.sim$n.sims / covid.sim$n.chains)))
  }
  
  params_names <- grepl("deviance", names(covid.sim$sims.list))
  params_names <- names(covid.sim$sims.list)[!params_names]
  
  for (name in params_names){
    
    if(length(covid.sim$sims.list[name][[1]]) < (covid.sim$n.sims + 1)){
      
      print(paste0("Saving traceplot for ", name))
      col_name <- paste0(name)
      
      sims <- covid.sim$sims.list[name][[1]]
      sim <- sims
      
      p <- as_tibble(sim) |>
        mutate(chain = factor(chain), n_sim = n_sim) |>
        ggplot() +
        geom_line(aes(n_sim, value, color = chain)) +
        labs(y = col_name)
      
      png(paste0("../results/", model_id, "/", col_name, "_tp.png"))
      print(p)
      dev.off()
      
    }
    
    else{
     
      for (k in 1:dim(covid.sim$sims.list[name][[1]])[2]){
        
        print(paste0("Saving traceplot for ", name,"_",k))
        col_name <- paste0(name,"_",k)
        
        sims <- covid.sim$sims.list[name][[1]]
        sim <- sims[,k]
        
        p <- as_tibble(sim) |>
          mutate(chain = factor(chain), n_sim = n_sim) |>
          ggplot() +
          geom_line(aes(n_sim, value, color = chain)) +
          labs(y = col_name)
        
        png(paste0("../results/", model_id, "/", col_name, "_tp.png"))
        print(p)
        dev.off()
        
      }
       
    }
      
  }

}

save_mcmc_params <- function(model_id, covid.sim){
  
  
  mcmc_params <- list("n.chains" = covid.sim$n.chains, "n.iter" = covid.sim$n.iter,
                     "n.burnin" = covid.sim$n.burnin, "n.thin" = covid.sim$n.thin,
                     "n.keep" = covid.sim$n.keep, "n.sims" = covid.sim$n.sims)
  
  sink(paste0("../results/", model_id, "/mcmc_params.txt"))
  print(mcmc_params)
  sink()
  
}


############## Procesamiento inicial de los datos ##############

### Leemos los datos, filtramos por region, calculamos las 
### metricas de interes, y realizamos cualquier proceso
### adicional segun la region.

# Carga de datos
datos <-  read.csv("../data/covid.csv")

# Filtramos la region de interes
region <- "South America"
datos  <- datos |> filter(continent == region)

# Variable auxiliar para desfazar las defunciones
aux <- datos |>
  select(c(Fatalities)) |>
  lag() |>
  pull(Fatalities)

# Obtenemos defunciones esperadas
# !!! Aqui el 5% es un supuesto muy importante !!!
covid <- datos |>
  mutate(total_deaths = Population * Deathrate / 1000000) |>
  mutate(expected_covid_deaths = (5 / 100 * total_deaths))

# Calculamos defunciones diarias observadas
# !!! Aqui lo del -1 puede que genere malos datos si
# el ultimo del pais anterior fue menor al primero del pais
# que le sigue
covid <- covid |>
  mutate(Fatalities_per_day = covid |> pull(Fatalities) - aux) |>
  filter(Fatalities_per_day > -1)

# Calculamos las medias de las metricas
covid <- covid |>
  select(c("Country", "Fatalities_per_day", 
           "expected_covid_deaths", "Net.migration")) |>
  group_by(Country) |>
  summarise(mean_fatalities = mean(Fatalities_per_day),
            mean_deathrate = mean(expected_covid_deaths),
            mean_net_migration = mean(Net.migration)) |>
  mutate(mean_fatalities = round(mean_fatalities, 0) |> as.integer())

# Procesamiento adicional segun la region
if (region == "South America"){
  covid <- covid |> filter(Country != "Suriname")
}


######
# covid <- covid |> filter(Country != "Brazil")
# covid <- covid |> filter(Country != "Bolivia")
######

# Conjunto de datos final  
covid

######################### Coordenadas #########################

# Definimos la lista en la que guardaremos los poligonos
lista_poligonos <- list()
lista_coords_x <- list()
lista_coords_y <- list()

# Guardamos los poligonos en la lista que definimos
for (pais in covid |> pull(Country)){
  lista_poligonos[pais] <- obtener_poligono(pais)
  lista_coords_x[pais] <- (lista_poligonos[pais][[1]]@Polygons[[1]]@coords |> colMeans())[1]
  lista_coords_y[pais] <- (lista_poligonos[pais][[1]]@Polygons[[1]]@coords |> colMeans())[2]
}

# Y creamos el objeto SpatialPolygons
sps = SpatialPolygons(lista_poligonos)
plot(sps)

# Guardamos las coordenadas centrales de cada pais
mean_coords <- tibble(
  Country = covid$Country,
  lista_coords_x |> as_tibble() |> t() |> as_tibble() |>  rename(x = V1), 
  lista_coords_y |> as_tibble() |> t() |> as_tibble() |>  rename(y = V1)
)

# Pequenios ajustes para que no se encimen
mean_coords <- mean_coords |>
  mutate(x = replace(x, Country == "Argentina", -64)) |>
  mutate(y = replace(y, Country == "Argentina", -32)) |>
  mutate(x = replace(x, Country == "Brazil", -52)) |>
  mutate(x = replace(x, Country == "Chile", -78)) |>
  mutate(y = replace(y, Country == "Chile", -39)) |>  
  mutate(x = replace(x, Country == "Colombia", -74)) |>
  mutate(x = replace(x, Country == "Ecuador", -86)) |>
  mutate(y = replace(y, Country == "Ecuador", 4)) |>
  mutate(x = replace(x, Country == "Guyana", -54)) |>
  mutate(y = replace(y, Country == "Guyana", 9)) |>
  mutate(x = replace(x, Country == "Paraguay", -78)) |>
  mutate(y = replace(y, Country == "Paraguay", -21)) |>
  mutate(x = replace(x, Country == "Peru", -76)) |>
  mutate(x = replace(x, Country == "Uruguay", -50)) |>
  mutate(y = replace(y, Country == "Uruguay", -36)) |>
  mutate(y = replace(y, Country == "Venezuela", 9)) |>
  mutate(x = replace(x, Country == "Venezuela", -66))
  
  
################# Graficas de las coordenadas #################

# Procesamiento adicional segun la region
if (region == "South America"){
  x_etiquetas <- -40
  y_etiquetas <- -15
  nclr <- 4
}

plot_arrows <- function(region){
  
  if (region == region){
    arrows(-75, -39, -72, -38, length = 0.1) #Chile
    arrows(-56, 8, -59, 5, length = 0.1) #Guyana
    arrows(-81, 2, -78, -1, length = 0.1) #Ecuador
    arrows(-75, -22.5, -59, -22.5, length = 0.1) #Paraguay
    arrows(-52, -35, -56, -32, length = 0.1) #Uruguay
  }
  
}

  #### SMR (standarized morbility rate)
# En este caso, analizamos lambda = observados / esperados
{
  plotvar <- covid$mean_fatalities/covid$mean_deathrate
  plotvar <- ifelse(plotvar == Inf, 0, plotvar)
  plotclr <- brewer.pal(nclr,"YlOrBr")
  class <- classIntervals(plotvar, nclr, dataPrecision=2,style="quantile")
  colcode <- findColours(class,plotclr)
  
  covid.map <- sps
  
  png(paste0("../results/smr.png"))
  plot(covid.map, col = colcode)
  legend(x_etiquetas, y_etiquetas, legend = names(attr(colcode, "table")), 
         fill = attr(colcode, "palette"), cex=1, bty="n")
  text(mean_coords$x, mean_coords$y, mean_coords$Country)
  text(mean_coords$x, mean_coords$y-3, plotvar |> round(2))
  plot_arrows(region)
  title(main="SMR")
  dev.off()
}

#### SMR (standarized morbility rate)
# En este caso, analizamos lambda = log(observados) / log(esperados)
{
  plotvar <- log(covid$mean_fatalities+1)/log(covid$mean_deathrate+1)
  plotvar <- ifelse(plotvar == Inf, 0, plotvar)
  plotclr <- brewer.pal(nclr,"YlOrBr")
  class <- classIntervals(plotvar, nclr, dataPrecision=2,style="quantile")
  colcode <- findColours(class,plotclr)
  
  covid.map <- sps
  
  png(paste0("../results/smr_log_log.png"))
  plot(covid.map, col = colcode)
  legend(x_etiquetas, y_etiquetas, legend = names(attr(colcode, "table")), 
         fill = attr(colcode, "palette"), cex=1, bty="n")
  text(mean_coords$x, mean_coords$y, mean_coords$Country)
  text(mean_coords$x, mean_coords$y-3, plotvar |> round(2))
  plot_arrows(region)
  title(main="SMR")
  dev.off()
}

#### Net Migration
{
  plotvar <- covid$mean_net_migration
  plotclr <- brewer.pal(nclr,"YlOrBr")
  class <- classIntervals(plotvar, nclr, dataPrecision=2,style="quantile")
  colcode <- findColours(class,plotclr)
  
  covid.map <- sps

  png(paste0("../results/net_migration.png"))
  plot(covid.map, col = colcode)
  legend(x_etiquetas, y_etiquetas, legend = names(attr(colcode, "table")), 
         fill = attr(colcode, "palette"), cex=1, bty="n")
  text(mean_coords$x, mean_coords$y, mean_coords$Country)
  text(mean_coords$x, mean_coords$y-3, plotvar |> round(2))
  plot_arrows(region)
  title(main="Net Migration")
  dev.off()
}



################# Matrices de adyacencia #################

n <- nrow(covid)

W.nb <- poly2nb(sps)
print(W.nb)
W.nb_list <- head(W.nb, 100)
max_n_neighbours <- 0
for (pais in W.nb_list){
  if (max_n_neighbours < length(pais)){
    max_n_neighbours <- length(pais)
  }
}
m <- rep(0,n)
W.l <- matrix(NA, nrow=n , ncol=max_n_neighbours )
adj <- NULL

for (i in 1:n){
  if (W.nb[[i]][1]!=0){
    m[i] <- length(W.nb[[i]])
    W.l[i,1:m[i]] <- W.nb[[i]]
    adj <- c(adj,W.nb[[i]])
  }
}

W <- matrix(0, nrow = n, ncol = n)

for (i in 1:n){
  for (j in 1:m[i]){
    W[i,W.l[i,j]] <- 1
    W[W.l[i,j],i] <- 1
  }
}

weights <- rep(1,length(adj))

################### Parametros CAR Propio ###################

########## Verificacion de eigenvalores ########## 
#
# El calculo del los eigenvalores para definir el rango de
# rho difiere entre las notas de clase y el manual de OpenBugs.
# Se observo que esto se debe a que las matrices D (en Bugs es M) y
# W (en Bugs es C) difieren en la definicion.
#
# Esta subseccion tiene el proposito de verificar que ambas formas de
# calcular los eigenvalores sea la misma.
#
# Nota: Para calcular ambas se requiere que NO haya regiones aisladas.
#
# https://www.mrc-bsu.cam.ac.uk/wp-content/uploads/geobugs12manual.pdf
# https://www.multibugs.org/documentation/latest/spatial/SpatialDistributions.html#Conditional

#### Notas de Clase ####

D_neg_sqr <- diag(m^(-1/2))

B <- D_neg_sqr %*% W %*% D_neg_sqr
eigenvalues_clase <- eigen(B)$values
max(eigenvalues_clase)
min(eigenvalues_clase)

#### OpenBugsManual ####

D_sqr <- diag((1/m)^(1/2))
D_neg_sqr <- diag((1/m)^(-1/2))

W_OB <- matrix(W, nrow = n)
for (k in 1:n){
  W_OB[k,] <- W_OB[k,]/m[k]
}

B <- D_neg_sqr %*% W_OB %*% D_sqr
# B <- D_neg_sqr %*% W_OB %*% D_neg_sqr

eigenvalues_bugs <- eigen(B)$values
max(eigenvalues_bugs)
min(eigenvalues_bugs)

#### Comparacion ####

# La suma debe ser cercana a 0
sum((eigenvalues_clase - eigenvalues_bugs)^2)

########## Definicion de parametros ########## 
#
# Una vez verificado que los eigenvalores sean los mismos,
# ocupamos la version de Bugs para que corra.
#

mu_areas <- rep(0, n)

M <- 1/m
C <- matrix(W_OB)
C <- C[C != 0]

### rho puede estar en el rango (1/min_eigen, 1/max_eigen)
max_eigen <- max(eigenvalues_bugs)
min_eigen <- min(eigenvalues_bugs)
rho_weight <- 90 # Elegir del 0 al 100
rho <- (rho_weight/100) * (1 / max_eigen) + ((100-rho_weight)/100) * (1 / min_eigen)

############ Asignacion de variables para modelo ############
#
# !!! Importante: hay que elegir solo alguno de los siguientes 
# !!! chunks para definir y, ee
#

{
  y <- covid$mean_fatalities
  ee <- covid$mean_deathrate  
}

{
  y <- covid$mean_fatalities/1000
  ee <- covid$mean_deathrate/1000
}

{
  y <- log(covid$mean_fatalities + 1)
  ee <- log(covid$mean_deathrate + 1)
}

x <- covid$mean_net_migration

x <- (covid$mean_net_migration |> scale(T,T))[,1]

#-Defining data-

data_normal <- list("n" = n, "y" = y, "ee" = ee, "x" = x, 
             "adj" = adj, "weights" = weights, "num" = m)

data_proper <- list("n" = n, "y" = y, "ee" = ee, "x" = x, 
             "adj" = adj, "num" = m, "mu_areas" = mu_areas,
             "C"=C, "M"=M, "gamma"=rho)

data_proper_nox <- list("n" = n, "y" = y, "ee" = ee, 
                    "adj" = adj, "num" = m, "mu_areas" = mu_areas,
                    "C"=C, "M"=M, "gamma"=rho)

data_proper_rho_prior <- c("n" = n, "y" = y, "ee" = ee, "x" = x, 
                           "adj" = adj, "num" = m, "mu_areas" = mu_areas,
                           "C"=C, "M"=M,  #quitamos gamma
                           "min_ie"=(1 / min_eigen),"max_ie"=(1 / max_eigen))
# data_proper_rho_prior <- c(data_proper[-length(data_proper)]) #quitamos gamma

data_hierarchical <- list("n" = n, "y" = y, "ee" = ee, "x" = x)

#-Defining inits-

inits <- function(){list(beta = rep(0,2),
                       tau.t = 1,
                       tau.c = 1,
                       theta = rep(0,n),
                       phi = rep(0,n),
                       yf = rep(0,n))}

inits_proper_nox <- function(){list(beta = 0,
                         tau.t = 1,
                         tau.c = 1,
                         theta = rep(0,n),
                         phi = rep(0,n),
                         yf = rep(0,n))}

inits_proper_rho_prior <- function(){list(beta = rep(0,2),
                                         tau.t = 1,
                                         tau.c = 1,
                                         gamma.c = 0.85,
                                         theta = rep(0,n),
                                         phi = rep(0,n),
                                         yf = rep(0,n))}

inits_normal_normal <- function(){list(beta = rep(0,2),
                         tau.y = 1,
                         tau.t = 1,
                         tau.c = 1,
                         theta = rep(0,n),
                         phi = rep(0,n),
                         yf = rep(0,n))}

inits_hierarchical <- function(){list(beta = rep(0,2),
                         mu.t = 0, 
                         tau.t = 1,
                         theta = rep(0,n),
                         yf = rep(0,n))}

#-Selecting parameters to monitor-

parameters <- c("beta", "lambda", "theta", "phi", "yf")

parameters_hierarchical <- c("beta", "lambda", "theta", "yf")

#-Running code-

covid.sim_normal <- bugs(data_normal, inits, parameters, model.file = "covid_car_normal.txt",
                         n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                         debug = T)

covid.sim_proper <- bugs(data_proper, inits, parameters, model.file = "covid_car_proper.txt",
                         n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                         debug = T)

covid.sim_proper_nox <- bugs(data_proper_nox, inits_proper_nox, parameters, model.file = "covid_car_proper_nox.txt",
                         n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                         debug = T)
# Este no corre
covid.sim_proper_rho_prior <- bugs(data_proper_rho_prior, inits_proper_rho_prior, parameters,
                                   model.file = "covid_car_proper_rho_prior.txt",
                                   n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                                   debug = T)

covid.sim_normal_normal <- bugs(data_normal, inits_normal_normal, parameters,
                                model.file = "covid_car_normal_normal.txt",
                                n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                                debug = T)

covid.sim_proper_normal <- bugs(data_proper, inits_normal_normal, parameters,
                                model.file = "covid_car_proper_normal.txt",
                                n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                                debug = T)

covid.sim_hierarchical <- bugs(data_hierarchical, inits_hierarchical, parameters_hierarchical,
                                model.file = "covid_hierarchical.txt",
                                n.iter = 10000, n.chains = 3, n.burnin = 1000, n.thin = 1,
                                debug = T)

########### Identificacion del modelo ###########
#
# Elegimos un identificador para guardar la corrida actual

model_id <- "car_proper_rho_90_nox"

###### !!! Importante !!! ######
#
# Hay que correr solo 1 de los siguientes chunks, dependiendo
# de cual es el model que se acaba de correr.
#
# Tambien de preferenci hay que elegir un nombre de identificador
# diferente si es que cambiaron varias cosas.

### Modelo CAR Intrinseco
{
  covid.sim <- covid.sim_normal
  save_definition(model_id, data_normal, inits, parameters, covid.sim)
}

### Modelo CAR Propio
{
  covid.sim <- covid.sim_proper
  save_definition(model_id, data_proper, inits, parameters, covid.sim)
}

### Modelo CAR Propio (Sin x)
{
  covid.sim <- covid.sim_proper_nox
  save_definition(model_id, data_proper_nox, inits_proper_nox, parameters, covid.sim)
}

### Modelo CAR Intrinseco con log(y) Normal
{
  covid.sim <- covid.sim_normal_normal
  save_definition(model_id, data_normal, inits_normal_normal, parameters, covid.sim)
}

### Modelo CAR Propio con log(y) Normal
{
  covid.sim <- covid.sim_proper_normal
  save_definition(model_id, data_proper, inits_normal_normal, parameters, covid.sim)
}

### Modelo Hierarchical
{
  covid.sim <- covid.sim_hierarchical
  save_definition(model_id, data_hierarchical, inits_hierarchical,
                  parameters_hierarchical, covid.sim)
}

############ Verificacion del modelo ############

#Traza de la cadena
# traceplot(covid.sim)

#Cadena y resumen
{
  out <- covid.sim$sims.list
  out.sum <- covid.sim$summary
}

#DIC
{
  out.dic<-covid.sim$DIC
  print(out.dic)
}

#Beta
{
  z<-out$beta[,2]
  png(paste0("../results/", model_id, "/beta_2_summary.png"))
  par(mfrow=c(2,2))
  plot(z,type="l")
  plot(cumsum(z)/(1:length(z)),type="l")
  hist(z,freq=FALSE)
  acf(z)
  dev.off()
}


#Tabla resumen
{
  out.b <- out.sum[grep("beta",rownames(out.sum)),c(1,3,7)]
  out.b <- cbind(out.b,apply(out$beta,2,prob))
  dimnames(out.b)[[2]][4] <- "prob"
  print(out.b)
  sink(paste0("../results/", model_id, "/beta_summary.txt"))
  print(out.b)
  sink()
}

#Predictions
{
  out.yf <- out.sum[grep("yf",rownames(out.sum)),]
  or <- order(y)
  ymin <- min(y, out.yf[,c(1,3,7)])
  ymax <- max(y, out.yf[,c(1,3,7)])
}

{
  png(paste0("../results/", model_id, "/predictions.png"))
  par(mfrow=c(1,1))
  plot(y[or],ylim=c(ymin,ymax), xlab="", xaxt='n')
  axis(1, at=1:n, labels=covid$Country[or], las=2)
  lines(out.yf[or,1],lwd=2,col=2)
  lines(out.yf[or,3],lty=2,col=2)
  lines(out.yf[or,7],lty=2,col=2)
  dev.off()
}

#R^2
{
  png(paste0("../results/", model_id, "/r2.png"))
  par(mfrow=c(1,1))
  plot(y,out.yf[,1], ylab="Predicted y")
  text(y,out.yf[,1], covid$Country)
  abline(a=0, b=1)
  dev.off()
  R2<-(cor(y,out.yf[,1]))^2
  print(R2)
  sink(paste0("../results/", model_id, "/r2.txt"))
  print(R2)
  sink()
}

#phi
{
  out.phi<-out.sum[grep("phi",rownames(out.sum)),]
  out.est<-out.phi
  k<-n
  ymin<-min(out.est[,c(1,3,7)])
  ymax<-max(out.est[,c(1,3,7)])
}

{
  png(paste0("../results/", model_id, "/phi.png"))
  par(mfrow=c(1,1))
  plot(1:k,out.est[,1][or], ylab="",ylim=c(ymin,ymax),
       xlab="", xaxt='n')
  axis(1, at=1:n, labels=covid$Country[or], las=2)
  segments(1:k,out.est[,3][or],1:k,out.est[,7][or])
  abline(h=0,col="grey70")
  title("Spatial Effect")
  dev.off()
}

#theta
{
  out.the<-out.sum[grep("the",rownames(out.sum)),]
  out.est<-out.the
  k<-n
  ymin<-min(out.est[,c(1,3,7)])
  ymax<-max(out.est[,c(1,3,7)])
}

{
  png(paste0("../results/", model_id, "/theta.png"))
  par(mfrow=c(1,1))
  plot(1:k,out.est[,1][or],ylab="",ylim=c(ymin,ymax),
       xlab="", xaxt='n')
  axis(1, at=1:n, labels=covid$Country[or], las=2)
  segments(1:k,out.est[,3][or],1:k,out.est[,7][or])
  abline(h=0,col="grey70")
  title("Individual Effect")
  dev.off()
}


#Map of lambda
{
  out.lam<-out.sum[grep("lam",rownames(out.sum)),]
  plotvar<-out.lam[,1]
  plotclr <- brewer.pal(nclr,"YlOrBr")
  class <- classIntervals(plotvar, nclr, dataPrecision=2,style="quantile")
  colcode <- findColours(class,plotclr)
  
  covid.map <- sps

  png(paste0("../results/", model_id, "/predicted_smr.png"), res = 70)
  plot(covid.map, col = colcode)
  legend(x_etiquetas, y_etiquetas, legend = names(attr(colcode, "table")), 
         fill = attr(colcode, "palette"), cex=1, bty="n")
  text(mean_coords$x, mean_coords$y, mean_coords$Country)
  text(mean_coords$x, mean_coords$y-3, plotvar |> round(2))
  plot_arrows(region)
  title(main="Smoothed SMR")
  dev.off()
}

#### SMR (standarized morbility rate)
#
# Esta solo es para comparar rapidamente con la que se genera en el
# chunk anterior

{
  plotvar <- y / ee
  plotvar <- ifelse(plotvar == Inf, 0, plotvar)
  plotclr <- brewer.pal(nclr,"YlOrBr")
  class <- classIntervals(plotvar, nclr, dataPrecision=2,style="quantile")
  colcode <- findColours(class,plotclr)
  
  covid.map <- sps

  plot(covid.map, col = colcode)
  
  legend(x_etiquetas, y_etiquetas, legend = names(attr(colcode, "table")), 
         fill = attr(colcode, "palette"), cex=1, bty="n")
  text(mean_coords$x, mean_coords$y, mean_coords$Country)
  text(mean_coords$x, mean_coords$y-3, plotvar |> round(2))
  plot_arrows(region)
  title(main="SMR")
}

############## Resumen de resultados ##############

########## log(y+1),  log(e+1) ##########

log_log_folders <- list.files("../results/", "log_log",
                            recursive=TRUE, include.dirs=TRUE)

log_log_folders <- log_log_folders[!(log_log_folders |> grepl(pattern = ".png"))]

log_log_dics <- list()

for (folder in log_log_folders){
  
  dic <- read.table(paste0("../results/",folder,"/dic.txt"))
  log_log_dics[folder] <- dic[2][[1]]
  
}

log_log_dics_df <-  log_log_dics |> 
  as_tibble() |> 
  t() |>
  as_tibble() |>
  rename(dic = V1) |> 
  mutate(model  = log_log_folders |> 
           gsub(pattern = "_log_log", replacement = "")) |>
  select(c(model, dic)) |>
  arrange(dic)

log_log_dics_df

############## y, e ##############

folders <- list.files("../results/", "", include.dirs=TRUE)

folders <- folders[!(folders |> grepl(pattern = ".png"))]
folders <- folders[!(folders |> grepl(pattern = "_log_log"))]
folders <- folders[!(folders |> grepl(pattern = "_div"))]

dics <- list()

for (folder in folders){
  
  dic <- read.table(paste0("../results/",folder,"/dic.txt"))
  dics[folder] <- dic[2][[1]]
  
}

dics_df <-dics |> 
  as_tibble() |> 
  t() |>
  as_tibble() |>
  rename(dic = V1) |> 
  mutate(model  = folders)|>
  select(c(model, dic)) |>
  arrange(dic)

dics_df
