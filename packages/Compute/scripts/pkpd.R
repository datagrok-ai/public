#name: pkpd
#language: r
#input: double dosage = 1000
#input: double doseInterval = 12
#output: graphics PKPD
#tags: model

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

# Define system parameters
params <- 
  c(KA = 0.3, CL = 17,             # central 
    V2 = 4.0E+01, Q = 1.0E+01, V3 = 3E+02, # peripheral
    Kin = .2, Kout = .2, EC50 = 8)              # effects
inits <- c(0, 0, 0, 1)   

# Plot results

ev <- eventTable()

# Specify 5 days BID dosing
ev$add.dosing(dose = dosage, nbr.doses = floor(120/doseInterval), dosing.interval = doseInterval)

# Use “start.time” parameter to specify 5 days QD starting at the end of the 5 days BID period


# Simulate and plot
Sims <- as.data.frame(mod1$run(params, ev, inits))

cMax <- max(Sims$centr)
effMax <- max(Sims$eff)

ggplot() + 
  geom_line(data = Sims, aes(x = time, y = centr), col = "Red", alpha = 0.7) + theme_bw()+
  geom_line(data = Sims, aes(x = time, y = eff*0.5*cMax/effMax), col = "Blue", alpha = 0.7) + theme_bw() +
  scale_y_continuous(name = "Drug concentration, nM", 
                     sec.axis = sec_axis( trans=~.*effMax*2/cMax, name="Drug effect")) +
  scale_x_continuous(name = "Time, hours", breaks = seq(0, 120, 12))
  