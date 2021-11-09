#name: pkpd
#language: r
#input: double dosage = 1000 {category: Dosing options} 
#input: double doseInterval = 12 {category: Dosing options} 
#input: string compartments {category: PK model; choices: ['1 compartment PK', '2 compartment PK']} 
#input: double clearance = 17 {category: PK parameters} 
#input: double rateConstant = 0.3 {category: PK parameters} [rate constant] 
#input: double centralV = 4 {category: PK parameters} [central compartment volume]  
#input: double perV = 30 {category: PK parameters} [peripheral compartment volume]
#input: double interRate = 1 {category: PK parameters} [intercompartmental rate] 
#input: double effRate = 0.2 {category: PD parameters} [effective compartment rate]
#input: double effect = 8 {category: PD parameters} [EC50]
#output: graphics PKPD
#tags: model

require("ggplot2")
require("RxODE")
rxCreateCache()

params <- 
  c(KA = rateConstant, CL = clearance,          # central 
    V2 = centralV, Q = interRate, V3 = perV,      # peripheral
    Kin = effRate, Kout = effRate, EC50 = effect)              # effects
inits <- c(0, 0, 0, 1) 

ev <- eventTable()
ev$add.dosing(dose = dosage, nbr.doses = floor(120/doseInterval), dosing.interval = doseInterval)

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

Sims <- as.data.frame(mod$run(params, ev, inits))

cMax <- max(Sims$centr)
effMax <- max(Sims$eff)

ggplot() + 
  geom_line(data = Sims, aes(x = time, y = centr), col = "Red", alpha = 0.7) + theme_bw()+
  geom_line(data = Sims, aes(x = time, y = eff*0.5*cMax/effMax), col = "Blue", alpha = 0.7) + theme_bw() +
  scale_y_continuous(name = "Drug concentration, nM", 
                     sec.axis = sec_axis( trans=~.*effMax*2/cMax, name="Drug effect")) +
  scale_x_continuous(name = "Time, hours", breaks = seq(0, 120, 12))
  