library(formatR)
library(graph)
library(gridExtra)
library(INLA)
library(lattice)
library(latticeExtra)
library(maptools)
library(R2WinBUGS)
library(RColorBrewer)
library(Rgraphviz)
library(rgdal)
library(spdep)
library(viridis)

dir.create("r")
dir.create("figuras")
dir.create("informes")
dir.create("datos")
dir.create("datos/brutos")
dir.create("datos/procesados")

load(file.path("..", "datos", "procesados", "Aragon.Rdata"))
aragon.shp <- readOGR(file.path("..", "datos", "procesados", "aragon.shp"))

aragon.shp <- aragon.shp[order(aragon.shp$CODMUNI), ]
aragon.nb <- poly2nb(aragon.shp)
vecinos <- nb2WB(aragon.nb)
paleta <- colorRampPalette(brewer.pal(9, "Blues"))(5)

modelo <- function() {
    for (i in 1:n) {
        O[i] ~ dpois(mu[i])
        log(mu[i]) <- log(E[i]) + m + het[i] + sp[i]
        het[i] ~ dnorm(0, prechet)
        R[i] <- exp(m + het[i] + sp[i])
    }
    sp[1:n] ~ car.normal(adj[], w[], num[], precsp)
    m ~ dflat()
    prechet <- pow(sdhet, -2)
    precsp <- pow(sdsp, -2)
    sdhet ~ dunif(0, 10)
    sdsp ~ dunif(0, 10)
}

set.seed(123)

datos <- list(E = Aragon.df$E, O = Aragon.df$O, n = dim(Aragon.df)[1], adj = vecinos$adj, w = vecinos$weights, num = vecinos$num)
parametros <- c("mu", "sdhet", "sdsp", "m", "R")
iniciales <- function() {
    list(sdhet = runif(1), sdsp = runif(1), m = rnorm(1))
}
iteraciones <- 20000
burnin <- 2000

ajuste.modelo <- bugs(model = modelo, inits = iniciales, data = datos, parameters.to.save = parametros, n.iter = iteraciones,
    n.burnin = burnin)
save(ajuste.modelo, file = "ajuste.modelo.rda")
load(file.path("..", "r", "ajuste.modelo.rda"))



plot(aragon.shp, col = paleta[findInterval(ajuste.modelo$mean$R, c(0, 0.6, 0.9, 1, 1.1, 1.8))])
title("SMR ")
legend("bottomright", c("0 - 0.6", "0.6 - 0.9", "0.9 - 1", "1 - 1.1", "1.1 - 1.8"), fill = paleta, cex = 0.7)



Nareas <- length(Aragon.df[, 1])
temp <- poly2nb(aragon.shp)
nb2INLA(file.path("..", "figuras", "LDN.graph"), temp)

H <- inla.read.graph(filename = file.path("..", "figuras", "LDN.graph"))
S <- U <- seq(1, 729)
data <- cbind(Aragon.df, S, U)
formula <- O ~ 1 + f(S, model = "besag", graph = H, scale.model = TRUE, hyper = list(prec = list(prior = "loggamma",
    param = c(1, 0.001)))) + f(U, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.001))))
modelo.inla <- inla(formula, family = "poisson", data = Aragon.df, E = E, control.compute = list(dic = TRUE, waic = TRUE,
    cpo = TRUE), control.predictor = list(compute = TRUE, cdf = c(log(1))))
summary(modelo.inla)

aragon.shp$SMR_mean <- modelo.inla$summary.fitted.values$mean  # Media
aragon.shp$SMR_sd <- modelo.inla$summary.fitted.values$sd  # Desviación tipica
aragon.shp$SMR_p1 <- 1 - modelo.inla$summary.fitted.values$`1 cdf`  # Probabilidad de ser mayor que 1

SMR.cutoff <- c(0.6, 0.9, 1, 1.1, 1.8)
SMR_p1.cutoff <- c(0, 0.2, 0.8, 1)
SMR_disc <- cut(aragon.shp$SMR_mean, breaks = SMR.cutoff, include.lowest = TRUE)
SMR_p1_disc <- cut(aragon.shp$SMR_p1, breaks = SMR_p1.cutoff, include.lowest = TRUE)

aragon.shp$SMR_disc <- SMR_disc
aragon.shp$SMR_p1_disc <- SMR_p1_disc
grid.arrange(spplot(aragon.shp, c("SMR_disc"), col.regions = brewer.pal(9, "Blues")[c(2, 4, 6, 8)], main = "SMR", par.settings = list(axis.line = list(col = "transparent"))),
    spplot(aragon.shp, c("SMR_p1_disc"), col.regions = brewer.pal(9, "Blues")[c(3, 6, 9)], main = "p(SMR > 1)", par.settings = list(axis.line = list(col = "transparent"))),
    ncol = 2)
