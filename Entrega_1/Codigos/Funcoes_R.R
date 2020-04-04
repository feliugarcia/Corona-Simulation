if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "MASS", "plyr","data.table", "gganimate", "Rfast", "gifski", "ggplot2", "grid", "png", "maps")


atualizar.direcao <- function(dados) {
  # Definindo a direção e a nova posição dos que estão em casa
  ind.casa <- which(dados$casa==T)
  ind.sair <- ind.casa[sample(c(T, F), length(ind.casa), 
                              prob=c(p.sair.casa, 1-p.sair.casa),
                              replace=T)]
  # Calcula o angulo (em radianos) entre a posição atual e o centro da cidade
  angulo.aux <- atan2(centro[1] - dados$lat[ind.sair], 
                      centro[2] - dados$long[ind.sair])
  dados$direc[ind.sair] <- rnorm(length(ind.sair), angulo.aux, angulo.sigma)
  
  
  # Definindo a situação dos que estão fora de casa em movimento
  ind.fora.move <- which(dados$casa==F & dados$direc!=0)
  aux.fora.move <- ind.fora.move[sample(c(T, F), 
                                        length(ind.fora.move), 
                                        prob=c(p.parar.fora, 1-p.parar.fora),
                                        replace=T)]
  dados$direc[aux.fora.move] <- 0  # Parando os individuos selecionados
  
  
  # Definindo a situação dos que estão parados fora de casa
  ind.fora.parado <- setdiff(which(dados$casa==F), ind.fora.move)
  aux.fora.parado <- ind.fora.parado[sample(c(T, F), 
                                            length(ind.fora.parado), 
                                            prob=c(p.voltar.casa, 1-p.voltar.casa),
                                            replace=T)]
  # Calcula o angulo (em radianos) entre a posição atual e o centro da cidade
  angulo.aux <- atan2(casas$lat[aux.fora.parado] - dados$lat[aux.fora.parado], 
                      casas$long[aux.fora.parado] - dados$long[aux.fora.parado])
  dados$direc[aux.fora.parado] <- angulo.aux
  
  moveram <- c(ind.sair, setdiff(ind.fora.move, aux.fora.move), aux.fora.parado)
  moveram.saudaveis <- intersect(moveram, dados[dados$status==0, "ind"])
  i<-1
  aux <- (0:(n.grid-1))*n.grid
  for(ind in moveram.saudaveis) {
    indices.lat.grid <- which.min(abs(dados[ind, "lat"] - grid.lat))
    indices.lat <- aux + indices.lat.grid
    indice.long <- which.min(abs(dados[ind, "long"] - superficie[indices.lat, 2]))
    linha.grid[ind] <<- indices.lat[indice.long]
    i <- i+1
  }
  
  return(dados$direc)
}


atualizar.posicao <- function(dados) {
  ### Atualizar a posição dos indivíduos
  delta <- tamanho.passo * cbind(sin(dados$direc), 
                                 cos(dados$direc))
  delta[dados$direc==0, ] <- 0  # casos em que a direção é nula não devem se mover 
  return(dados[, c("lat", "long")] + delta)
}


atualizar.casa <- function(dados) {
  ### Atualizar a coluna casa
  distancia.casa <- sqrt(rowSums((dados[, c("lat", "long")] - casas[, c("lat", "long")])^2))
  estao.casa <- which(distancia.casa < tamanho.passo)
  dados[estao.casa, c("lat", "long")] <- casas[estao.casa, c("lat", "long")]
  dados$casa <- (1:n.pop) %in% estao.casa
  return(dados$casa)
}


atualizar.risco <- function(dados) {
  ### Atualizar a superficie de risco
  tmp <- dados %>%
    filter(status==1)
  aux <- kde2d(tmp$long, tmp$lat, h=rep(raio.risco, 2), n=n.grid, 
               lims=c(range(long.extremos), range(lat.extremos)))
  risco <- superficie[, tempo+1]*taxa.decaimento + as.vector(aux$z)
  return(risco)
}


atualizar.status <- function(dados) {
  ### Atualizar a coluna status
  # Simulando infecções
  risco <- superficie[, tempo+2]  # Nível de contaminação do ambiente
  p.infec <- 2*((1/(1+exp(-escala.risco*risco)))-.5)  # Probabilidade de infecção
  
  saudaveis <- dados[dados$status==0, "ind"]
  infectados <- as.integer((runif(length(saudaveis)) < p.infec[linha.grid[saudaveis]]))
  dados$status[saudaveis][infectados] <- 1L

  # Simulando recuperações
  contaminados <- setdiff(1:n.pop, saudaveis)
  aux <- sample(c(T, F), length(contaminados), 
                    prob=c(p.recuperar, 1-p.recuperar),
                    replace=T)
  recuperados <- contaminados[aux]
  dados$status[recuperados] <- 2L
  return(dados$status)
}
