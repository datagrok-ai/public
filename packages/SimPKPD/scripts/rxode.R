#name: rxodeCommandReal
#language: r
#input: double inputSD = 1
#output: graphics plot

require("ggplot2")
require("RxODE")
rxCreateCache()


mod1 <-RxODE({
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
})

params <- 
   c(KA = 0.3, CL = 7,             # central 
     V2 = 4.0E+01, Q = 1.0E+01, V3 = 3E+02, # peripheral
     Kin = .2, Kout = .2, EC50 = 8)              # effects
inits <- c(0, 0, 0, 1)   

# Initialize event table
ev <- eventTable()
# Specify dose
ev$add.dosing(dose = 10000, nbr.doses = 1)
# Specify sampling
ev$add.sampling(0:240)
# Simulate
Sims <- data.frame(mod1$run(params, ev, inits))
# Plot results
ggplot() + geom_line(data = Sims, aes(x = time, y = centr), col = "Red", alpha = 0.7) + theme_bw()