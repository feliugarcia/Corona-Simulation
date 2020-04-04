### Carregando os pacotes necessários
source("Funcoes_R.R")
set.seed(12345) # Definindo a semente

### Definindo os parâmetros do modelo
TT <- 12*24*7 # Número de passos (tempo 5 min)

### Parametros da cidade
n.pop <- 4600 # Tamanho da população (divisível por 8)
centro <- c(lat=-20.902805, long=-45.127847)
lat.extremos <- c(lat.min=-20.908992, lat.max=-20.897918)
long.extremos <- c(long.min=-45.132165, long.max=-45.123509)

### Parametros de deslocamento
p.sair.casa <- .05 # Probabilidade de uma pessoa em casa sair no minuto em questão
p.voltar.casa <- .05 # Probabilidade de uma pessoa parada fora de casa entrar em movimento de retorno
# A probabilidade abaixo induz uma distribuição geométrica para o "número de movimentos até parar"
p.parar.fora <- .2 # Probabilidade de uma pessoa em movimento fora de casa parar
# Desvio padrão do ruído que desloca a direção de movimentação em relação ao centro da cidade
angulo.sigma <- 1  # Quanto maior o valor menor é a chance dos moradores irem para o centro
# .0001 é aproximadamente 30 metros
tamanho.passo <- .0002 # Velocidade do deslocamento por unidade de tempo 


### Parametros de disseminação do vírus
taxa.decaimento <- .99  # Taxa com que a contaminação do ambiente é preservada no tempo seguinte
raio.risco <- .0005  # parâmetro que define o risco de contaminação ao redor de um indivíduo
escala.risco <- 1E-8  # valor que relaciona o risco (nível de contaminação) à probabilidade de infecção 
n.grid <- 100  # Número de pontos no grid de estimação da superficie de risco
p.recuperar <- .00001  # Probabilidade de uma pessoa infectada se recuperar em um instante de tempo


### Criando a lista que irá receber os dados simulados
dados <- vector(TT, mode="list")
names(dados) <- 1:(TT) 

### Definindo os locais de residência
casas <- data.frame(
  ind=1:n.pop,
  lat=runif(n.pop, lat.extremos[1], lat.extremos[2]),
  long=runif(n.pop, long.extremos[1], long.extremos[2])
)
# Induzir pessoas a morar juntas
casas[1:n.pop, ] <- casas[sample(n.pop/3, n.pop, replace=T), ]  


### Criando a população no tempo 1
dados[[1]] <- data.frame(
  ind=1:n.pop,
  lat=c(casas$lat[1:(n.pop/2)], 
        runif(n.pop/2, lat.extremos[1], lat.extremos[2])),
  long=c(casas$long[1:(n.pop/2)], 
         runif(n.pop/2, long.extremos[1], long.extremos[2])),
  casa=c(rep(1L, n.pop/2), rep(0L, n.pop/2)),
  status=c(rep(1L, 1), rep(0L, n.pop-1)), # 0=saudável, 1=doente, 2=recuperado
  # trocar a linha abaixo para mover as pessoas como na função atualizar.direcao 
  direc=rep(0, n.pop) 
) 
indices <- (n.pop/8*7):n.pop 
angulo.aux <- atan2(centro[1] - dados[[1]]$lat[indices], 
                    centro[2] - dados[[1]]$long[indices])
dados[[1]]$direc[indices] <- rnorm(length(indices), angulo.aux, angulo.sigma)


### Criando a superficie de risco
superficie <- matrix(nr=n.grid^2, nc=2+TT) # Superficie de risco
colnames(superficie) <- c("long", "lat", paste0("densidade_", 1:TT))
tmp <- dados[[1]] %>%
  filter(status==1)
aux <- kde2d(tmp$long, tmp$lat, h=rep(raio.risco, 2), n=n.grid, 
             lims=c(range(long.extremos), range(lat.extremos)))
superficie[, 1:2] <- unlist(with(aux, expand.grid(x, y)))
superficie[, 3] <- as.vector(aux$z)
grid.lat <- unique(superficie[, 1])
grid.long <- unique(superficie[, 2])


### Definindo o índice da linha correspondente ao ponto do grid mais próximo à posição corrente
linha.grid <- numeric(n.pop)
for(ind in 1:n.pop) {
  indices.lat.grid <- which.min(abs(dados[[1]][ind, "lat"] - grid.lat))
  indices.lat <- (0:(n.grid-1))*n.grid + indices.lat.grid
  indice.long <- which.min(abs(dados[[1]][ind, "long"] - superficie[indices.lat, 2]))
  linha.grid[ind] <- indices.lat[indice.long]
}


### Simulando a evolução da população
for(tempo in 2:TT) {
  dados[[tempo]] <- dados[[tempo-1]]
  dados[[tempo]]$direc <- atualizar.direcao(dados[[tempo]])
  dados[[tempo]][, c("lat", "long")] <- atualizar.posicao(dados[[tempo]])
  dados[[tempo]]$casa <- atualizar.casa(dados[[tempo]])
  superficie[, 2+tempo] <- atualizar.risco(dados[[tempo]])
  dados[[tempo]]$status <- atualizar.status(dados[[tempo]])
}


### Sumarizando os resultados
dados.plot <- as.data.frame(data.table::rbindlist(dados, idcol="tempo"))
class(dados.plot$tempo) <- "integer"
dados.plot %>%
  group_by(tempo) %>%
  dplyr::summarise(saudaveis=sum(status==0),
                   infectados=sum(status==1),
                   recuperados=sum(status==2)) %>%
  melt(id.vars="tempo", variable.name="status") %>%
  ggplot(aes(tempo, value, color=status)) +
  geom_line()

save.image("Dados.Rdata")

### Animação
anim.dinamica <- dados.plot %>% 
  filter(tempo<=288*2, ind %in% c(1, sample(n.pop, 199))) %>% 
  ggplot(aes(long, lat)) + 
  # geom_point(aes(color=as.factor(status), shape=as.factor(casa), size=status + .1)) +
  geom_point(aes(color=as.factor(status), shape=as.factor(casa))) +
  # facet_wrap(~ status) +
  transition_time(tempo) +
  ease_aes('cubic-in-out')
anim_save("dinamica.gif", animate(anim.dinamica, nframes=288*2, fps=10))



### Superfície
anim.superficie <- superficie %>%
  as_tibble() %>%
  melt(id.vars=c("lat", "long"), value.name="risco") %>%
  mutate(tempo=as.integer(rep(1:TT, each=nrow(superficie)))) %>%
  filter(tempo<=288*2) %>%
  ggplot(aes(long, lat, alpha=risco)) +
  geom_raster() +
  transition_time(tempo) +
  ease_aes('linear')
anim_save("superficie.gif", animate(anim.superficie, nframes=288*2, fps=10))
