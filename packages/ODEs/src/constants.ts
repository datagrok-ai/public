// Control constants
export const CONTROL_TAG = '#';
export const DF_NAME = 'df';

/** Control expressions for the problem specifying */
export enum CONTROL_EXPR {
    NAME = `${CONTROL_TAG}name`,
    TAGS = `${CONTROL_TAG}tags`,
    DESCR = `${CONTROL_TAG}description`,
    DIF_EQ = `${CONTROL_TAG}equations`,
    EXPR = `${CONTROL_TAG}expressions`,
    ARG = `${CONTROL_TAG}argument`,
    INITS = `${CONTROL_TAG}inits`,
    CONSTS = `${CONTROL_TAG}constants`,
    PARAMS = `${CONTROL_TAG}parameters`,
    TOL = `${CONTROL_TAG}tolerance`,
    LOOP = `${CONTROL_TAG}loop`,
    UPDATE = `${CONTROL_TAG}update`,
};

/** Loop consts */
export enum LOOP {
  MIN_LINES_COUNT = 1,
  COUNT_IDX = 0,
  COUNT_NAME = '_count',  
  MIN_COUNT = 1,
};

/** UPDATE consts */
export enum UPDATE {
  MIN_LINES_COUNT = 1,
  DURATION_IDX = 0,
  DURATION = '_duration',  
  MIN_DURATION = 0,
};

/** Basic template illustrating the simplest features */
const TEMPLATE_BASIC = `${CONTROL_EXPR.NAME}: Template 
${CONTROL_EXPR.DIF_EQ}:
  dy/dt = -y + sin(t) / t

${CONTROL_EXPR.ARG}: t
  initial = 0.01
  final = 15
  step = 0.01

${CONTROL_EXPR.INITS}:  
  y = 0`;

/** Advanced template illustrating extanded features */
const TEMPLATE_ADVANCED = `NOTES. This is an advanced template. Modify it. Use multi-line formulas if needed.
Add new equations, expressions, constants & parameters. Edit these notes if required.

${CONTROL_EXPR.NAME}: Advanced
${CONTROL_EXPR.DIF_EQ}:
  dx/dt = E1 * y + sin(t)

  dy/dt = E2 * x - pow(t, 5)

${CONTROL_EXPR.EXPR}:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2

${CONTROL_EXPR.ARG}: t
  start = 0
  finish = 2
  step = 0.01

${CONTROL_EXPR.INITS}:
  x = 2
  y = 0

${CONTROL_EXPR.CONSTS}:
  C1 = 1
  C2 = 3

${CONTROL_EXPR.PARAMS}:
  P1 = 1
  P2 = -1

${CONTROL_EXPR.TOL}: 0.00005`;

/** Extended template illustrating extanded features */
const TEMPLATE_EXTENDED = `NOTES. This is an extended template. It has additional scripting annotations.
Use it as a backbone for models & apps creating.

${CONTROL_EXPR.NAME}: Extended
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: 2D ordinary differential equations system sample
${CONTROL_EXPR.DIF_EQ}:
  dx/dt = E1 * y + sin(t)
  dy/dt = E2 * x - pow(t, 5)

${CONTROL_EXPR.EXPR}:
  E1 = C1 * exp(-t) + P1
  E2 = C2 * cos(2 * t) + P2

${CONTROL_EXPR.CONSTS}:
  C1 = 1
  C2 = 3

${CONTROL_EXPR.PARAMS}:
  P1 = 1 {category: Parameters} [P1 parameter]
  P2 = -1 {category: Parameters} [P2 parameter]

${CONTROL_EXPR.INITS}:
  x = 2 {units: C; category: Initial values} [Initial value of x]
  y = 0 {units: C; category: Initial values} [Initial value of y]

${CONTROL_EXPR.ARG}: t
  start = 0 {units: min; caption: Initial; category: Time} [Initial time of simulation]
  finish = 2 {units: min; caption: Final; category: Time} [Final time of simulation]
  step = 0.01 {units: min; caption: Initial; category: Time} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.00005`;

/** Chemical reactions, mass-action kinetics */
const CHEM_REACT_MODEL = `The following example is taken from https://doi.org/10.1002/ijch.201800003.

${CONTROL_EXPR.NAME}: Chem react
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Mass-action kinetics illustration
${CONTROL_EXPR.DIF_EQ}:
  dx1/dt = -k1 * x1 + k2 * (x2)**2 + k3 * x1 * x3 
           - k4 * (x1)**2 - 2 * k5 * (x1)**2 + k6 * x2 * x4

  dx2/dt = 2 * k1 * x1 - 2 * k2 * (x2)**2 
           + k5 * (x1)**2 - k6 * x2 * x4

  dx3/dt = -k3 * x1 * x3 + k4 * (x1)**2 + k6 * x2 * x4

  dx4/dt = k5 * (x1)**2 - k6 * x2 * x4

${CONTROL_EXPR.PARAMS}:
  k1 = 0.7 {category: Reaction parameters}
  k2 = 0.9 {category: Reaction parameters}
  k3 = 1.2 {category: Reaction parameters}
  k4 = 3.4 {category: Reaction parameters}
  k5 = 2.3 {category: Reaction parameters}
  k6 = 4.5 {category: Reaction parameters}

${CONTROL_EXPR.INITS}:
  x1 = 1 {units: mol/L; category: Initial concentrations}
  x2 = 0 {units: mol/L; category: Initial concentrations}
  x3 = 0 {units: mol/L; category: Initial concentrations}
  x4 = 0 {units: mol/L; category: Initial concentrations}

${CONTROL_EXPR.ARG}: t
  initial = 0 {units: min; caption: Initial; category: Time} [Initial time of simulation]
  final = 5 {units: min; caption: Final; category: Time} [Final time of simulation]
  step = 0.01 {units: min; caption: Step; category: Time} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.00005`;

/** Robertson's chemical reaction model - stiff ODEs */
const ROBERTSON_MODEL = `NOTES. The classic example of stiff ODEs: the Robertson problem.

${CONTROL_EXPR.NAME}: Robertson
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Robertson's chemical reaction model
${CONTROL_EXPR.DIF_EQ}:
  dA/dt = -0.04 * A + 10000 * B * C
  dB/dt = 0.04 * A - 10000 * B * C - 30000000 * B**2
  dC/dt = 30000000 * B**2

${CONTROL_EXPR.INITS}:
  A = 1 {units: mol/L; category: Initial concentrations}
  B = 0 {units: mol/L; category: Initial concentrations}
  C = 0 {units: mol/L; category: Initial concentrations}

${CONTROL_EXPR.ARG}: t
  start = 0 {units: sec; caption: Initial; category: Time} [Initial time of simulation]
  finish = 40 {units: sec; caption: Final; category: Time} [Final time of simulation]
  step = 0.01 {units: sec; caption: Step; category: Time} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.0000001`;

/** Fermentation process simulation */
const FERMENTATION_MODEL = `NOTES. The following problem is taken from https://core.ac.uk/download/pdf/11737483.pdf.

${CONTROL_EXPR.NAME}: Fermentation
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Simulation of fermentation process in the ethanol production
${CONTROL_EXPR.DIF_EQ}:
  dP/dt = r * X
  dS/dt = -q * X
  dX/dt = V * S / (K + S) * X

${CONTROL_EXPR.ARG}: t
  initial = 0 {units: d; caption: Initial; category: Time} [Initial time of simulation]
  final = 70 {units: d; caption: Final; category: Time} [Final time of simulation]
  step = 0.01 {units: d; caption: Step; category: Time} [Time step of simlulation]

${CONTROL_EXPR.INITS}:  
  P = 4.276 {units: ml/ml; caption: ethanol (P); category: Initials} [Concentration of ethanol]
  S = 0.3185 {units: mg/ml; caption: glucose (S); category: Initials} [Concentration of glucose]
  X = 0.092 {units: mg/ml; caption: Saccharomyces (X); category: Initials} [Saccharomyces wet weight]

${CONTROL_EXPR.PARAMS}:
  r = 1.3455 {units: mg/ml; category: Parameters} [The growth rate of ethanol]
  q = 0.011129 {units: mg/ml; category: Parameters} [The rate of glucose consumption]
  V = 0.087 {units: mg/ml; category: Parameters} [The maximal growth rate of Saccharomyces]
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
  ${LOOP.COUNT_NAME} = 10 {category: Dosing} [Number of doses]
  depot += dose

${CONTROL_EXPR.ARG}: t
  start = 0 {units: h; caption: begin; category: Dosing} [Begin of dosing interval]
  final = 12 {units: h; caption: end; category: Dosing} [End of dosing interval]
  step = 0.01 {units: h; caption: step; category: Dosing} [Time step of simlulation]  

${CONTROL_EXPR.INITS}:  
  depot = 0 {category: Initial values}
  centr = 0 {category: Initial values} [Central]
  peri = 0 {category: Initial values} [Peripheral]
  eff = 0.2 {category: Initial values} [Effective compartment rate]

${CONTROL_EXPR.PARAMS}:  
  dose = 10000 {category: Dosing} [Dosage]
  KA = 0.3 {caption: rate constant; category: Paramters}
  CL = 2 {caption: clearance; category: Paramters}
  V2 = 4 {caption: central volume; category: Paramters} [Central compartment volume]
  Q = 1 {caption: intercompartmental rate; category: Paramters}
  V3 = 30 {caption: peripheral volume; category: Paramters} [Peripheral compartment volume]
  EC50 = 8 {caption: effect; category: Paramters}
  Kin = 0.2 {caption: Kin; category: Paramters} [The first-order production constant]
  Kout = 0.2 {caption: Kout; category: Paramters} [The first-order dissipation rate constant]
  
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
  stage1 = 60 {units: h; caption: 1-st stage; category: Time} [Duration of the 1-st stage]
  step = 0.1 {units: h; caption: step; category: Time} [Time step of simlulation]

${CONTROL_EXPR.UPDATE}:
  stage2 = 60 {units: h; caption: 2-nd stage; category: Time} [Duration of the 2-nd stage]
  S += 70

${CONTROL_EXPR.INITS}:  
  X = 5 {units: kg/m^3; caption: biomass; category: Initial concentrations} [Aspergillus niger biomass]
  S = 150 {units: kg/m^3; caption: glucose; category: Initial concentrations} [Glucose]
  O = 7 {units: kg/m^3; caption: oxygen; category: Initial concentrations} [Dissolved oxygen]
  P = 0 {units: kg/m^3; caption: acid; category: Initial concentrations} [Gluconic acid]

${CONTROL_EXPR.PARAMS}:  
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

/** Initial value problem templates */
export enum TEMPLATES {
  EMPTY = '',
  BASIC = TEMPLATE_BASIC,
  ADVANCED = TEMPLATE_ADVANCED,
  EXTENDED = TEMPLATE_EXTENDED,
  CHEM_REACT = CHEM_REACT_MODEL,
  ROBERTSON = ROBERTSON_MODEL,
  FERMENTATION = FERMENTATION_MODEL,
  PK_PD = PK_PD_MODEL,
  ACID_PROD = ACID_PROD_MODEL
};