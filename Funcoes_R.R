if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse", "MASS", "plyr","data.table", "gganimate", "Rfast", "gifski", "readr", "tibble", "plotrix")


atualizar.direcao <- function(dados) {
  
  ### Atualizar a direção de movimentação.
  
  # Definindo a direção e a nova posição dos que estão em casa
  ind.casa <- which(dados$casa==TRUE)
  ind.sair <- ind.casa[sample(c(TRUE, FALSE), length(ind.casa), 
                              prob=c(p.sair.casa, 1-p.sair.casa),
                              replace=TRUE)]
  
  # Calcula o angulo (em radianos) entre a posição atual e o centro da cidade
  
  desloc_func <- sample(dim(locais)[[1]],length(ind.sair), replace = TRUE)
  angulo.aux <- atan2(locais[desloc_func,2] - dados$lat[ind.sair],
                      locais[desloc_func,3] - dados$long[ind.sair])
  
  dados$direc[ind.sair] <- rnorm(length(ind.sair), angulo.aux, angulo.sigma)
  
  # Definindo a situação dos que estão fora de casa em movimento
  ind.fora.move <- which(dados$casa==FALSE & dados$direc!=0)
  
  aux.fora.move <- ind.fora.move[sample(c(TRUE, FALSE), 
                                        length(ind.fora.move), 
                                        prob=c(p.parar.fora, 1-p.parar.fora),
                                        replace=TRUE)]
  dados$direc[aux.fora.move] <- 0  # Parando os individuos selecionados
  
  # Definindo a situação dos que estão parados fora de casa e não hospitalizados ou mortos
  ind.fora.parado <- setdiff(which(dados$casa==FALSE & 
                                     dados$status<=2 & 
                                     dados$resid==TRUE), 
                             ind.fora.move)
  
  aux.fora.parado <- ind.fora.parado[sample(c(TRUE, FALSE), 
                                            length(ind.fora.parado), 
                                            prob=c(p.voltar.casa, 1-p.voltar.casa),
                                            replace=TRUE)]
  
  
  # Calcula o angulo (em radianos) entre a posição atual e o centro da cidade
  angulo.aux <- atan2(casas$lat[aux.fora.parado] - dados$lat[aux.fora.parado], 
                      casas$long[aux.fora.parado] - dados$long[aux.fora.parado])
  
  dados$direc[aux.fora.parado] <- angulo.aux
  

  #Movimenta os saudáveis, inclusive os forasteiros.
  
  moveram <- c(ind.sair, setdiff(ind.fora.move, aux.fora.move), aux.fora.parado)
  moveram.saudaveis <- intersect(moveram, dados[dados$status==0, "ind"])
  i<-1
  aux <- (0:(n.grid-1))*n.grid
  for(ind in moveram.saudaveis) {
    if (ind > n.pop){
      
      amostra <- as.vector(dados[dados$ind<=n.pop & dados$status==0, "ind"])
      if (length(amostra) > 0){
        aux.ind <- sample(amostra, size = 1, replace = TRUE)
        indices.lat.grid <- which.min(abs(dados[aux.ind, "lat"] - grid.lat))
        indices.lat <- aux + indices.lat.grid
        indice.long <- which.min(abs(dados[aux.ind, "long"] - superficie[indices.lat, 2]))
        linha.grid[ind] <<- indices.lat[indice.long]
      }
      
    }else{
      
      indices.lat.grid <- which.min(abs(dados[ind, "lat"] - grid.lat))
      indices.lat <- aux + indices.lat.grid
      indice.long <- which.min(abs(dados[ind, "long"] - superficie[indices.lat, 2]))
      linha.grid[ind] <<- indices.lat[indice.long]
      
    }
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
  
  tam.dados <- length(dados[, "ind"]) # dados$resid==1
  
  ### Atualizar a coluna casa

  ## Atualiza as casas para os residentes.
  tmp <- dados[dados$resid==1, c("lat", "long")]
  d <- (tmp - casas[, c("lat", "long")])
  distancia.casa <- sqrt(rowSums(d)^2)
  estao.casa <- which(distancia.casa < tamanho.passo)
  dados[estao.casa, c("lat", "long")] <- casas[estao.casa, c("lat", "long")]
  dados$casa <- (1:tam.dados) %in% estao.casa
  
  ## Atualiza as casas dos forasteiros
  if (tam.dados > n.pop){
    
    tmp <- dados[dados$resid==0, c("lat", "long")]
    t <- length(tmp[, "lat"])
    casas.forast.index <- sample(c(1:length(casas[, "ind"])), t, replace = TRUE)
    d <- (tmp - c(casas$lat[casas.forast.index], casas$long[casas.forast.index]))
    distancia.casa <- sqrt(rowSums(d)^2)
    estao.casa <- which(distancia.casa < tamanho.passo)
    if (length(estao.casa) > 0){
      dados[estao.casa, c("lat", "long")] <- casas[estao.casa, c("lat", "long")]
      dados$casa[estao.casa] <- TRUE 
    }
  }
  
  
  return(dados$casa)
}


atualizar.risco <- function(dados, tempo, superficie) {
  
  ### Atualizar a superficie de risco
  tmp <- dados %>% filter(status==1)
  if (empty(tmp)){
    risco <- superficie[, tempo+1]*taxa.decaimento
  }else{
    aux <- kde2d(na.omit(tmp$long), na.omit(tmp$lat), h=rep(raio.risco, 2), n=n.grid, 
                 lims=c(range(long.extremos), range(lat.extremos)))
    risco <- superficie[, tempo+1]*taxa.decaimento + as.vector(aux$z)
  }
  return(risco)
}

#########################################################################################
# Função responsável atualização do status da população.
# Este método aplica as probabilidades de infecção, hospitalização, intenração em UTI,
# óbito e recuperação.
#########################################################################################
atualizar.status <- function(dados, tempo, superficie, linha.grid) {
  
  tam.dados <- length(dados[, "ind"])
  
  ### Atualizar a coluna status
  # Simulando infecções
  risco <- superficie[, tempo+2]  # Nível de contaminação do ambiente
  p.infec <- 2*((1/(1+exp(-escala.risco*risco)))-.5)  # Probabilidade de infecção

  p.infec.forast <- mean(p.infec)

  saudaveis <- dados[dados$status==0 & dados$resid==1, "ind"]
  infectados <- as.integer((runif(length(saudaveis)) < p.infec[linha.grid[saudaveis]]))
  dados$status[saudaveis][infectados] <- 1L

  saudaveis.forast <- dados[dados$status==0 & dados$resid==0, "ind"]
  aux_saudaveis_forast <- sample(c(TRUE, FALSE), length(saudaveis.forast), prob = c(p.infec.forast, 1-p.infec.forast), replace = TRUE)
  contaminados <- saudaveis.forast[aux_saudaveis_forast]
  if (length(contaminados) > 0){
    dados$status[dados$ind %in% contaminados] <- 1L 
  }

  # Silumando hospitalizações
  contaminados <- dados[dados$status==1, "ind"]
  aux_hospitalizados <- sample(c(TRUE, FALSE), length(contaminados),
                               prob=c(p.hospt, 1-p.hospt),
                               replace=TRUE)
  hospitalizados <- contaminados[aux_hospitalizados] 
  if (length(hospitalizados) > 0){
    dados$status[dados$ind %in% hospitalizados] <- 3L 
    dados$casa[dados$ind %in% hospitalizados] <- FALSE 
    dados$direc[dados$ind %in% hospitalizados] <- 0
  }  

  # Simulando intenrações em UTI
  hospitalizados <- dados[dados$status==3, "ind"]
  aux_uti <- sample(c(TRUE, FALSE), length(hospitalizados),
                    prob=c(p.uti, 1-p.uti),
                    replace=TRUE)
  uti <- hospitalizados[aux_uti] 
  if (length(uti) > 0){
    dados$status[dados$ind %in% uti] <- 4L 
    dados$casa[dados$ind %in% uti] <- FALSE 
    dados$direc[dados$ind %in% uti] <- 0
  }  
  
  # Simulando Óbitos
  uti <- dados[dados$status==4, "ind"]
  aux_obito <- sample(c(TRUE, FALSE), length(uti),
                      prob=c(p.obito, 1-p.obito),
                      replace=TRUE)
  obito <- uti[aux_obito] 
  if (length(obito) > 0){
    dados$status[dados$ind %in% obito] <- 5L 
    dados$casa[dados$ind %in% obito] <- FALSE 
    dados$direc[dados$ind %in% obito] <- 0
  }  

  # Simulando recuperações

  saudaveis <- dados[dados$status==0, "ind"]
  obitos <- dados[dados$status==5, "ind"]
  recuperados <- dados[dados$status==2, "ind"]

  contaminados <- setdiff(1:tam.dados, saudaveis)
  contaminados <- setdiff(contaminados, obitos)
  contaminados <- setdiff(contaminados, recuperados)

  aux <- sample(c(TRUE, FALSE), length(contaminados),
                prob=c(p.recuperar, 1-p.recuperar),
                replace=TRUE)
  recuperados <- contaminados[aux]
  dados$casa[recuperados] <- TRUE
  dados$status[recuperados] <- 2L
  
  return(dados$status)
}


#########################################################################################
# Função responsável pela adição de forasteiros à população.
# A seleção dos forasteiros utiliza um tamanho fixo, definido em configuração
# que é utilizado para definir um valor máximo em torno do qual o número de  forasteiros 
# variará ao longo do tempo.
#########################################################################################
adicionar.forasteiros <- function(dados, tempo, superficie){
  
  dados.interno <- dados
  
  ### Atualizar a população volante de forasteiros
  
    insereForasteiros <- sample(c(TRUE, FALSE), prob = c(get.p.forasteiros.entrar(), 1-get.p.forasteiros.entrar()))
  
    if (insereForasteiros && max.pop.forasteiros > 0){
    
    ind.forasteiros.existentes <- dados.interno[dados.interno$resid==0, "ind"]
    tam.pop.forasteiros <- 1:max.pop.forasteiros
    qtd.forasteiros.entra <- round(runif(1, 1, max.pop.forasteiros))
    id.forasteiros.entra <- (n.pop+1):(n.pop+qtd.forasteiros.entra)
    ind.forasteiros.entra <- setdiff(id.forasteiros.entra, ind.forasteiros.existentes)
    
    # Simulando infecções
    risco <- mean(sample(superficie[, tempo+2], qtd.forasteiros.entra))   # Nível de contaminação do ambiente
    p.infec <- 2*((1/(1+exp(-escala.risco*risco)))-.5)  # Probabilidade de infecção
    
    
    if (length(ind.forasteiros.entra) > 0)  {
      
      # criando o dataframe para incluir na população.
      forasteiros.entra <- data.frame(
        ind=ind.forasteiros.entra,
        lat=runif(length(ind.forasteiros.entra), lat.extremos[1], lat.extremos[2]),
        long=runif(length(ind.forasteiros.entra), long.extremos[1], long.extremos[2]),
        casa=0, # define que forasteiros não ocuparão casa.
        status=sample(x = c(1, 0),
                      size = length(ind.forasteiros.entra),
                      prob = c(p.infec, 1-p.infec),
                      replace = TRUE), # 0=saudável, 1=doente, 2=recuperado, 3=hospitalizado, 4=uti, 5=óbito
        # trocar a linha abaixo para mover as pessoas como na função atualizar.direcao
        direc=rep(0, length(ind.forasteiros.entra)),
        resid=rep(0, length(ind.forasteiros.entra)) # define que se trata de um forasteiro.
      )
      
      dados.interno <- rbind(dados.interno, forasteiros.entra)
      
    }
    
  }
  
  return (dados.interno[order(dados.interno$ind, decreasing = FALSE), ])
  
}


#########################################################################################
# Função responsável pela remoção de forasteiros à população.
# Este método faz com que a população de forasteiros varie ao longo do tempo.
#########################################################################################
remover.forasteiros <- function(dados){
  
  dados.interno <- dados
  
  if (max.pop.forasteiros > 0){

    excluirForasteiros <- sample(c(TRUE, FALSE), prob = c(get.p.forasteiros.sair(), 1-get.p.forasteiros.sair()))
    
    qtd.forasteiros.sai <- round(runif(1, 1, max.pop.forasteiros))
    
    ind.forasteiros.existentes <- dados.interno[dados.interno$resid==0, "ind"] 
    
    if (excluirForasteiros && length(ind.forasteiros.existentes) > qtd.forasteiros.sai){
      
      id.forasteiros.sai <- sample(ind.forasteiros.existentes, qtd.forasteiros.sai, replace = TRUE)
      ind.forasteiros.sai <- setdiff(ind.forasteiros.existentes, id.forasteiros.sai)
      
      dados.interno <- dados.interno[-ind.forasteiros.sai, ]

    }
  }
  
  return(dados.interno[order(dados.interno$ind, decreasing = FALSE), ])
}

### Simulando a evolução da população
executar_simulacao <- function(n.pop, TT){

  for(tempo in 2:TT) {
    
    #cat("Tempo: ", tempo, " - Início !\n")
    dados[[tempo]] <- dados[[tempo-1]]
    
    dados[[tempo]]$direc <- atualizar.direcao(dados[[tempo]])
    dados[[tempo]][, c("lat", "long")] <- atualizar.posicao(dados[[tempo]])

    dados[[tempo]]$casa <- atualizar.casa(dados[[tempo]])
    superficie[, 2+tempo] <- atualizar.risco(dados[[tempo]], tempo, superficie)

    dados[[tempo]] <- adicionar.forasteiros(dados[[tempo]], tempo, superficie)

    dados[[tempo]]$status <- atualizar.status(dados[[tempo]], tempo, superficie, linha.grid)

    dados[[tempo]] <- remover.forasteiros(dados[[tempo]])
    
  }
  cat("Tempo: ", tempo, " - Fim !\n")
  return(dados)
  
}
