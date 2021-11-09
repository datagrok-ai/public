#name: popPK
#language: r
#input: double clearance = 1
#input: double omegaClearance = 0.5
#output: graphics popPK
#tags: model

require("ggplot2")
require("RxODE")
require("dplyr")
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
  c(KA = 0.3, CL = clearance,             # central 
    V2 = 4.0E+01, Q = 1, V3 = 3E+02, # peripheral
    Kin = .2, Kout = .2, EC50 = 8)              # effects
inits <- c(0, 0, 0, 1)   

# Initialize event table
ev <- eventTable()
# Specify dose
ev$add.dosing(dose = 1000, nbr.doses = 1)

# Specify sampling
ev$add.sampling(0:240)

SimsAll <- data.frame()

for(i in seq(1, 20)){
  a <- params
  a["CL"] <- params["CL"] + rnorm(1, 0, omegaClearance)
  
  Sims <- data.frame(mod1$run(a, ev, inits)) %>% mutate(ID = i)
  SimsAll <- SimsAll %>% bind_rows(Sims)
}


DataSet <- SimsAll %>% select(time, centr, ID) %>% 
  group_by(time) %>% 
  summarise(mean = mean(centr), 
            q95 = quantile(centr, probs = 0.95),
            q05 = quantile(centr, probs = 0.05),
            n = n()) %>% 
  ungroup() 

# Plot results
ggplot() + 
  geom_line(data = DataSet, aes(x = time, y = mean), col = "#483D8B", size =1) + theme_bw()+
  geom_ribbon(data = DataSet, aes(x = time, ymin = q05, ymax = q95), col = "#ADFF2F", fill = "#ADFF2F", alpha = 0.5) + 
  scale_y_continuous(name = "Drug concentration, uM")+
  scale_x_continuous(name = "Time, hours", breaks = seq(0,1200, 24), limits = c(0, 240))