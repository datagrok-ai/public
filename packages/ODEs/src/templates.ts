import {CONTROL_EXPR} from './constants';

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
const TEMPLATE_ADVANCED = `${CONTROL_EXPR.NAME}: Advanced
${CONTROL_EXPR.COMMENT}:
  This is an advanced template. Modify it. Use multi-line formulas if needed.
  Add new equations, expressions, constants & parameters. Edit these comment lines if required.
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
const TEMPLATE_EXTENDED = `${CONTROL_EXPR.NAME}: Extended
${CONTROL_EXPR.TAGS}: model
${CONTROL_EXPR.DESCR}: 2D ordinary differential equations system sample
${CONTROL_EXPR.COMMENT}:
  This is an extended template. It has additional scripting annotations.
  Use it as a backbone for the platform applications creating.

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

/** Initial value problem templates */
export enum TEMPLATES {
  EMPTY = '',
  BASIC = TEMPLATE_BASIC,
  ADVANCED = TEMPLATE_ADVANCED,
  EXTENDED = TEMPLATE_EXTENDED,
};

/** Demo template */
export const DEMO_TEMPLATE = [`${CONTROL_EXPR.NAME}: My model`,
`${CONTROL_EXPR.DIF_EQ}:`,
`  dx/dt = y + sin(t)`,
`  dy/dt = x - exp(t)`,  
``,
`${CONTROL_EXPR.ARG}: t`,
`  initial = 0`,
`  final = 4`,
`  step = 0.01`,
``,
`${CONTROL_EXPR.INITS}:`,
`  x = 2`,
`  y = 1`];
