library(xtable)

set.seed(1101985)
n <- 50
tau <- 5

C <- round(runif(n, 0, 10), 1)
T <- round(runif(n, 0, 10), 1)
V <- pmin(C, T)
delta <- ifelse(T <= C, 1, 0)
delta <- delta[order(V)]
V <- V[order(V)]
minVtau <- pmin(V, rep(tau, n))
E <- ifelse(delta == 1 & V < tau, 1, 0)
E <- ifelse(delta == 0 & V < tau, NA, E)
KM.cens <- survfit(Surv(V, 1-delta) ~1)
surv.prob <- summary(KM.cens, times = minVtau)$surv
omega <- 1/surv.prob*ifelse(is.na(E)==TRUE, 0, 1)
X <- rbinom(n, 1, 0.5)
sum(ifelse(omega == 0, 0, omega*E*X))/sum(ifelse(omega == 0, 0, omega*X))
sum(ifelse(omega == 0, 0, omega*E*(1-X)))/sum(ifelse(omega == 0, 0, omega*(1-X)))
data.xtable <- data.frame(X, V, delta, E, minVtau, surv.prob, omega)


print(xtable(data.xtable, digits = c(0,0,1, 0, 0, 1, 2, 2)), include.rownames = FALSE )
