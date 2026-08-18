#R script for using ODE of the sequential equations in SDDM

library(deSolve)
system <- function(t, state, parms) {
x <- state[1]
y <- state[2]
dx <- y
dy <- -x
list(c(dx, dy))
}
times <- seq(0, 10, by = 0.01)

sol <- ode(
y = c(x = 1, y = 0),
times = times,
func = system,
parms = NULL
)

plot(sol[, "time"], sol[, "x"], type = "l")


#Dt Function for metformin
#Dt=function(x) {if (x > 50) return (1.0) else 0.0}  # treatment starts at t=50\n"

system_SDDM <- function(t,state,parms) {
with(as.list(c(state, parms)),{
G <- state[1]
I <- state[2]
B <- state[3]
V <- state[4]
S <- state[5]
A <- state[6]
Ic <- state[7]

dG=Ra+ Rhep * (1 - xiA * A)- Uii- k1 * (1 + etaA * A) * (I / RI) * G

secretion=B * sigma * (1 / (1 + rhoA * A)) * (G^n / (G^n + G0^n))
dI= -k2 * I + secretion

dB= B * rB * (G / Gnorm - 1) - dB * B + chiA * A * B

dIc= epsilon * (I - Ic)

glucose_excess = max(G - Gnorm, 0)
damage = alpha * glucose_excess + gamma * Ic
repair = delta * (1 - V) * (1 + zetaA * A)
dV= repair - damage * V
hazard = kappa * (1 - V) + theta * (G - Gnorm)^2
dS= hazard * (1 - S)

dA= kon * Dt * (1 - A) - koff * A   #metformin Dt is not used in dG state

list(c(dG,dI,dB,dV,dS,dA,dIc))
})

}

times <- seq(0, 365, by = 1)

init=y0=c(6.0,15.0,1.0,0.9,0.0,0.0,15.0)
pars=c(Ra=0.5,
       Rhep=0.7,
       Uii=0.9,
       k1=0.004,
       RI=2.0,
       # Insulin
       k2=0.3,
       sigma=8,
       G0=7.0,
       n=2.0,
   
       #Beta cells
       rB=0.0002,
       dB=0.0001,
       chiA=0.0005,
       # Chronic insulin
       epsilon=0.005,
    
       # Vascular
       zetaA=0.3,
       alpha=0.0005,
       gamma=0.00005,
       delta=0.005,
    
       # Stroke
       kappa=0.00005,
       theta=0.00002,
    
       # AMPK,
         kon=0.02,
         koff=0.01,
         Dt=0.5,
   
       #Modifiers
        etaA=0.4,
        xiA=0.4,
        rhoA=0.3,
        Gnorm=5.5)
   

sol <- ode(
y =y0,
times = times,
func = system_SDDM,
parms =pars,method="lsoda")

plot(sol,main=c("Glucose","Insulin","β-Cell Mass / Function","Vascular Health","Stroke Risk","AMPK Activation")








