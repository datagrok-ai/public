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

/** Use case 1: chemical reactions, mass-action kinetics */
const CHEM_REACT_MODEL = `The following example is taken from https://doi.org/10.1002/ijch.201800003.

${CONTROL_EXPR.NAME}: Chem react
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Illustration of mass-action kinetics
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

/** Use case 1: Robertson's chemical reaction model - stiff ODEs */
const ROBERTSON_MODEL = `NOTES. The classic example of stiff ODEs: the Robertson problem.

${CONTROL_EXPR.NAME}: Robertson
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: Robertson's chemical reaction model
${CONTROL_EXPR.DIF_EQ}:
  dx/dt = -0.04 * x + 10000 * y * z
  dy/dt = 0.04 * x - 10000 * y * z - 30000000 * y**2
  dz/dt = 30000000 * y**2

  ${CONTROL_EXPR.INITS}:
  x = 1 {units: mol/L; caption: X; category: Initial concentrations}
  y = 0 {units: mol/L; caption: Y; category: Initial concentrations}
  z = 0 {units: mol/L; caption: Z; category: Initial concentrations}

  ${CONTROL_EXPR.ARG}: t
  start = 0 {units: sec; caption: Initial; category: Time} [Initial time of simulation]
  finish = 40 {units: sec; caption: Final; category: Time} [Final time of simulation]
  step = 0.01 {units: sec; caption: Step; category: Time} [Time step of simlulation]

  ${CONTROL_EXPR.TOL}: 0.0000001`;

/** Case 3 */
const USE_CASE3 = 'NOTES. Super cool model.';

/** Initial value problem templates */
export enum TEMPLATES {
  EMPTY = '',
  BASIC = TEMPLATE_BASIC,
  ADVANCED = TEMPLATE_ADVANCED,
  EXTENDED = TEMPLATE_EXTENDED,
  CHEM_REACT = CHEM_REACT_MODEL,
  ROBERTSON = ROBERTSON_MODEL,
  CASE3 = USE_CASE3,
};