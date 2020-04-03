chaveminha <- 'AIzaSyBQ_nLiWG3Nyltmso4CCelVikkSKzIBMzc'



library(ggmap)


register_google(key=chaveminha,account_type = 'standard')

jacare <- c(lon = -45.130833333333, lat = -20.905)

santana <- get_map(jacare, zoom = 16, scale = 2)

numpessoas <- 500

pessoas <- data.frame(lon = rnorm(numpessoas, mean = -45.130833333, sd = .003), 
                     lat = rnorm(numpessoas, mean = -20.905, sd = .003), cor = sample(1:3,numpessoas, replace = TRUE) )
pessoas$cor <- as.factor(pessoas$cor)

ggmap(santana) + geom_point(aes(lon, lat, color = cor), data = pessoas)

from <- c(pessoas[1,1],pessoas[1,2])
to <- c(pessoas[10,1],pessoas[10,2])

rota <-route(from, to, mode = 'walking')

rota <-route(from, to, structure = 'route')

ggmap(santana) +
  geom_path(
    aes(x = lon, y = lat),  colour = "red",
    size = 1.5, alpha = .5,
    data = rota, lineend = "round"
  )
ls()






