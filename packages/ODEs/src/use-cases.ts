import {CONTROL_EXPR} from './constants';

/** Chemical reactions, mass-action kinetics */
const CHEM_REACT_MODEL = `${CONTROL_EXPR.NAME}: Chem react
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Mass-action kinetics illustration
${CONTROL_EXPR.COMMENT}: This model is taken from https://doi.org/10.1002/ijch.201800003.
${CONTROL_EXPR.DIF_EQ}:
  dx1/dt = -k1 * x1 + k2 * (x2)**2 + k3 * x1 * x3 
           - k4 * (x1)**2 - 2 * k5 * (x1)**2 + k6 * x2 * x4

  dx2/dt = 2 * k1 * x1 - 2 * k2 * (x2)**2 
           + k5 * (x1)**2 - k6 * x2 * x4

  dx3/dt = -k3 * x1 * x3 + k4 * (x1)**2 + k6 * x2 * x4

  dx4/dt = k5 * (x1)**2 - k6 * x2 * x4

${CONTROL_EXPR.PARAMS}:
  k1 = 0.7 {category: Reaction parameters; min: 0.1; max: 5}
  k2 = 0.9 {category: Reaction parameters; min: 0.1; max: 5}
  k3 = 1.2 {category: Reaction parameters; min: 0.1; max: 5}
  k4 = 3.4 {category: Reaction parameters; min: 0.1; max: 5}
  k5 = 2.3 {category: Reaction parameters; min: 0.1; max: 5}
  k6 = 4.5 {category: Reaction parameters; min: 0.1; max: 5}

${CONTROL_EXPR.INITS}:
  x1 = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x2 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x3 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}
  x4 = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 2}

${CONTROL_EXPR.ARG}: t
  initial = 0 {units: min; caption: Initial; category: Time; min: 0; max: 5} [Initial time of simulation]
  final = 6 {units: min; caption: Final; category: Time; min: 6; max: 10} [Final time of simulation]
  step = 0.01 {units: min; caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.00005`;

/** Robertson's chemical reaction model - stiff ODEs */
const ROBERTSON_MODEL = `${CONTROL_EXPR.NAME}: Robertson
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Robertson's chemical reaction model
${CONTROL_EXPR.COMMENT}: This is classic example of stiff ODEs.
${CONTROL_EXPR.DIF_EQ}:
  dA/dt = -0.04 * A + 10000 * B * C
  dB/dt = 0.04 * A - 10000 * B * C - 30000000 * B**2
  dC/dt = 30000000 * B**2

${CONTROL_EXPR.INITS}:
  A = 1 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  B = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}
  C = 0 {units: mol/L; category: Initial concentrations; min: 0; max: 5}

${CONTROL_EXPR.ARG}: t
  start = 0 {units: sec; caption: Initial; category: Time; min: 0; max: 1} [Initial time of simulation]
  finish = 40 {units: sec; caption: Final; category: Time; min: 2; max: 50} [Final time of simulation]
  step = 0.01 {units: sec; caption: Step; category: Time; min: 0.01; max: 0.1} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.0000001`;

/** Fermentation process simulation */
const FERMENTATION_MODEL = `${CONTROL_EXPR.NAME}: Fermentation
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Simulation of fermentation process in the ethanol production
${CONTROL_EXPR.COMMENT}: This problem is taken from https://core.ac.uk/download/pdf/11737483.pdf.
${CONTROL_EXPR.DIF_EQ}:
  dP/dt = r * X
  dS/dt = -q * X
  dX/dt = V * S / (K + S) * X

${CONTROL_EXPR.ARG}: t
  initial = 0 {units: d; caption: Initial; category: Time; min: 0; max: 10} [Initial time of simulation]
  final = 70 {units: d; caption: Final; category: Time; min: 15; max: 100} [Final time of simulation]
  step = 0.01 {units: d; caption: Step; category: Time; min: 0.01; max: 1} [Time step of simlulation]

${CONTROL_EXPR.INITS}:  
  P = 4.276 {units: ml/ml; caption: ethanol, P; category: Initials; min: 3; max: 6} [Concentration of ethanol]
  S = 0.3185 {units: mg/ml; caption: glucose, S; category: Initials; min: 0.1; max: 0.5} [Concentration of glucose]
  X = 0.092 {units: mg/ml; caption: saccharomyces, X; category: Initials; min: 0.01; max: 0.2} [Saccharomyces wet weight]

${CONTROL_EXPR.PARAMS}:
  r = 1.3455 {units: mg/ml; category: Parameters; min: 1; max: 2} [The growth rate of ethanol]
  q = 0.011129 {units: mg/ml; category: Parameters; min: 0.01; max: 0.2} [The rate of glucose consumption]
  V = 0.087 {units: mg/ml; category: Parameters; min: 0.01; max: 0.2} [The maximal growth rate of Saccharomyces]
  K = 0.0628 {category: Parameters} [Michaelis-Menten constant]
  
${CONTROL_EXPR.TOL}: 0.0000001`;

/** PK-PD simulation */
const PK_PD_MODEL = `${CONTROL_EXPR.NAME}: PK-PD
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Pharmacokinetic-pharmacodynamic (PK-PD) simulation: two-compartment model
${CONTROL_EXPR.DIF_EQ}:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
  d(peri)/dt  = Q * C2 - Q * C3
  d(eff)/dt  = Kin - Kout * (1 - C2/(EC50 + C2)) * eff

${CONTROL_EXPR.EXPR}:
  C2 = centr / V2
  C3 = peri / V3

${CONTROL_EXPR.LOOP}:
  count = 10 {category: Dosing; min: 1; max: 20} [Number of doses]
  depot += dose

${CONTROL_EXPR.ARG}: t
  start = 0 {units: h; caption: begin; category: Dosing; min: 0; max: 1} [Begin of dosing interval]
  final = 12 {units: h; caption: end; category: Dosing; min: 5; max: 15} [End of dosing interval]
  step = 0.01 {units: h; caption: step; category: Dosing; min: 0.01; max: 0.1} [Time step of simlulation]  

${CONTROL_EXPR.INITS}:  
  depot = 0 {category: Initial values}
  centr = 0 {category: Initial values} [Central]
  peri = 0 {category: Initial values} [Peripheral]
  eff = 0.2 {category: Initial values} [Effective compartment rate]

${CONTROL_EXPR.PARAMS}:  
  dose = 10000 {category: Dosing; min: 1000; max: 20000; step: 1000} [Dosage]
  KA = 0.3 {caption: rate constant; category: Paramters; min: 0.1; max: 1}
  CL = 2 {caption: clearance; category: Paramters; min: 1; max: 5}
  V2 = 4 {caption: central volume; category: Paramters; min: 1; max: 10} [Central compartment volume]
  Q = 1 {caption: inter rate; category: Paramters; min: 0.1; max: 1} [Intercompartmental rate]
  V3 = 30 {caption: peri volume; category: Paramters; min: 20; max: 40} [Peripheral compartment volume]
  EC50 = 8 {caption: effect; category: Paramters; min: 1; max: 10}
  Kin = 0.2 {caption: Kin; category: Paramters; min: 0.1; max: 0.5} [The first-order production constant]
  Kout = 0.2 {caption: Kout; category: Paramters; min: 0.1; max: 0.5} [The first-order dissipation rate constant]
  
${CONTROL_EXPR.TOL}: 0.000000001`;

/** Gluconic acid production */
const ACID_PROD_MODEL = `${CONTROL_EXPR.NAME}: GA-production
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Gluconic acid (GA) production by Aspergillus niger modeling
${CONTROL_EXPR.DIF_EQ}:
  dX/dt = rX
  dS/dt = -gamma * rX - lambda * X
  dO/dt = Kla * (Cod - O) - delta * rX - phi * X
  dP/dt = alpha * rX + beta * X

${CONTROL_EXPR.EXPR}:
  mu = muM * S / (Ks + S) * O / (Ko + O)
  rX = mu * X

${CONTROL_EXPR.ARG}: t
  start = 0 {units: h; caption: initial; category: Time} [Start of the process]
  stage1 = 60 {units: h; caption: 1-st stage; category: Time; min: 40; max: 80} [Duration of the 1-st stage]
  step = 0.1 {units: h; caption: step; category: Time; min: 0.01; max: 1} [Time step of simlulation]

${CONTROL_EXPR.UPDATE}: 2-nd stage
  duration = stage2
  S += 70

${CONTROL_EXPR.INITS}:  
  X = 5 {units: kg/m^3; caption: biomass; category: Initial concentrations; min: 1; max: 10} [Aspergillus niger biomass]
  S = 150 {units: kg/m^3; caption: glucose; category: Initial concentrations; min: 50; max: 200} [Glucose]
  O = 7 {units: kg/m^3; caption: oxygen; category: Initial concentrations; min: 1; max: 10} [Dissolved oxygen]
  P = 0 {units: kg/m^3; caption: acid; category: Initial concentrations; min: 0; max: 0.1} [Gluconic acid]

${CONTROL_EXPR.PARAMS}:
  stage2 = 60 {units: h; caption: 2-nd stage; category: Time; min: 40; max: 80} [Duration of the 2-nd stage]
  muM = 0.668 {units: 1/h; category: Parameters} [Monod type model parameter]
  alpha = 2.92 {category: Parameters} [Monod type model parameter]
  beta = 0.131 {units: 1/h; category: Parameters} [Monod type model parameter]
  gamma = 2.12 {category: Parameters} [Monod type model parameter]
  lambda = 0.232 {units: 1/h; category: Parameters} [Monod type model parameter]
  delta = 0.278 {category: Parameters} [Monod type model parameter]
  phi = 0.00487 {units: 1/h; category: Parameters} [Monod type model parameter]
  Ks = 130.9 {units: g/L; category: Parameters} [Monod type model parameter]
  Ko = 0.000363 {units: g/L; category: Parameters} [Monod type model parameter]
  Kla = 0.017 {units: 1/s; category: Parameters} [Volumetric mass transfer coefficient]
  Cod = 15 {units: kg/m^3; category: Parameters} [Liquid phase dissolved oxygen saturation concentration]
  
${CONTROL_EXPR.TOL}: 0.000000001`;

/** Nimotuzumab disposition model */
const NIMOTUZUMAB_MODEL = `${CONTROL_EXPR.NAME}: Nimotuzumab
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Nimotuzumab disposition model
${CONTROL_EXPR.COMMENT}: Source: https://www.mdpi.com/1999-4923/12/12/1147
${CONTROL_EXPR.DIF_EQ}:
  dA1/dt = (-(CL * A3 / V1 + Q / V1) * A1 + Q / V2 * A2 - kint * Rtot * A1 / (Kss + A1 / V1)) 
           / (1 + Rtot * Kss / (Kss + A1 / V1)**2)

  dA2/dt = (Q / V1 * A1 - Q / V2 * A2) / (1 + Rtotp * Kss / (Kss + A2 / V2)**2)

  dA3/dt = Kin * (1 + Smax * A1**gamma / (S50**gamma + A1**gamma)) - Kout * A3

${CONTROL_EXPR.ARG}: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10} [Initial time of simulation]
  finish = 70 {caption: Final; category: Time; min: 50; max: 100} [Final time of simulation]
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 1} [Time step of simlulation]

${CONTROL_EXPR.INITS}:
  A1 = 0.43 {category: Initials; min: 0.2; max: 0.6} [Central]
  A2 = 0.55 {category: Initials; min: 0.3; max: 0.8} [Periferal]
  A3 = 10.32 {category: Initials; min: 5; max: 15} [Mediator]

${CONTROL_EXPR.OUTPUT}:
  t {caption: Time, h}
  A1 {caption: Central}
  A2 {caption: Periferal}

${CONTROL_EXPR.PARAMS}:
  CL = 0.00964 {category: Parameters; min: 0.005; max: 0.015} [Non-specific clearance]
  V1 = 2.63 {category: Parameters; min: 0.7; max: 3} [Apparent volume of distribution of the central compartment]
  Q = 0.0288 {category: Parameters; min: 0.015; max: 0.035} [Distribution clearance of free nimotuzumab between the central and peripheral compartment]
  V2 = 0.00992 {category: Parameters; min: 0.0067; max: 0.0265} [Apparent volume of distribution of the peripheral compartment]
  Kss = 15.5 {category: Parameters; min: 10; max: 45} [Quasi steady state constant for the EGFR]
  kint = 0.01 {category: Parameters; min: 0.01; max: 0.05} [Internalization rate for nimotuzumab-EGFR]
  Rtot = 0.0105 {category: Parameters; min: 0.007; max: 0.032} [EGFR in the central compartment]
  Rtotp = 956 {category: Parameters; min: 120; max: 3500} [EGFR in the peripheral compartment]
  Kin = 0.0133 {category: Parameters; min: 0.01; max: 0.05} [Zero-order synthesis rate constant]
  Kout = 0.0133 {category: Parameters; min: 0.001; max: 0.02} [First-order elimination rate constant]
  S50 =	8.57 {category: Parameters; min: 5; max: 12} [Concentration of free nimotuzumab in the central compartment that achieves the half of Smax]
  Smax = 3.18 {category: Parameters; min: 1; max: 5} [Maximal effect of the stimulation]
  gamma = 0.5 {category: Parameters; min: 0.1; max: 1} [Hill coefficient of the sigmoid function]

${CONTROL_EXPR.TOL}: 0.0000001`;

/** Initial value problem use cases */
export enum USE_CASES {
  CHEM_REACT = CHEM_REACT_MODEL,
  ROBERTSON = ROBERTSON_MODEL,
  FERMENTATION = FERMENTATION_MODEL,
  PK_PD = PK_PD_MODEL,
  ACID_PROD = ACID_PROD_MODEL,
  NIMOTUZUMAB = NIMOTUZUMAB_MODEL,
};