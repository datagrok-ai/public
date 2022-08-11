#name: popPK
#language: r
#input: double dosage = 1000 {category: Dosing options} 
#input: double doseInterval = 12 {category: Dosing options} 
#input: string compartments {category: PK model; choices: ['2 compartment PK', '1 compartment PK']} 
#input: double clearance = 2 {category: PK parameters} 
#input: double omegaClearance = 0.3 {category: Random effects} [omega clearance] 
#input: double rateConstant = 0.3 {category: PK parameters} [rate constant] 
#input: double centralV = 4 {category: PK parameters} [central compartment volume]  
#input: double perV = 30 {category: PK parameters} [peripheral compartment volume]
#input: double interRate = 1 {category: PK parameters} [intercompartmental rate] 
#output: graphics popPK
#meta.domain: PKPD

require("ggplot2")
require("RxODE")
require("dplyr")
rxCreateCache()

params <- 
  c(KA = rateConstant, CL = clearance,          # central 
    V2 = centralV, Q = interRate, V3 = perV,      # peripheral
    Kin = effRate, Kout = effRate, EC50 = effect)              # effects
inits <- c(0, 0, 0, 1) 

# Initialize event table
ev <- eventTable()
# Specify dose
ev$add.dosing(dose = 1000, nbr.doses = 1)

if(compartments == "1 compartment PK"){
  mod <-RxODE({
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2;
    d/dt(peri)  = Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
  })
} else {
  mod <-RxODE({
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
  })
}


SimsAll <- data.frame()

for(i in seq(1, 20)){
  a <- params
  a["CL"] <- params["CL"] + rlnorm(1, 0, omegaClearance)
  
  Sims <- data.frame(mod$run(a, ev, inits)) %>% mutate(ID = i)
  SimsAll <- SimsAll %>% bind_rows(Sims)
}


DataSet <- SimsAll %>% select(time, centr, ID)

DataSet <- DataSet %>% 
  group_by(time) %>% 
  summarise(mean = mean(centr), 
            q95 = quantile(centr, probs = 0.95),
            q05 = quantile(centr, probs = 0.05),
            n = n()) %>% ungroup() 

# Plot results
ggplot() + 
  geom_line(data = DataSet, aes(x = time, y = mean), col = "#483D8B", size =1) + theme_bw()+
  geom_ribbon(data = DataSet, aes(x = time, ymin = q05, ymax = q95), col = "#ADFF2F", fill = "#ADFF2F", alpha = 0.5) + 
  scale_y_continuous(name = "Drug concentration, uM")+
  scale_x_continuous(name = "Time, hours", breaks = seq(0, 4, 24), limits = c(0, 24))+
  theme(axis.text.x = element_text(size = 24, angle = 0,vjust =.5, family="serif"),
        axis.text.y = element_text(size = 24, family = "serif"),
        axis.title.x = element_text(vjust = 1, size = 24, family = "serif"),
        axis.title.y = element_text(vjust = 1, size = 24, family = "serif")) 

omega <- omegaClearance