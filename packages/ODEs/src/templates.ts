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
  finish = 10
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
  P1 = 1 {category: Parameters; min: 1; max: 10} [P1 parameter]
  P2 = -1 {category: Parameters; min: -10; max: -1} [P2 parameter]

${CONTROL_EXPR.INITS}:
  x = 2 {category: Initial values; min: 0; max: 5} [Initial value of x]
  y = 0 {category: Initial values; min: -2; max: 2} [Initial value of y]

${CONTROL_EXPR.ARG}: t
  start = 0 {caption: Initial; category: Time; min: 0; max: 10} [Initial time of simulation]
  finish = 10 {caption: Final; category: Time; min: 10; max: 20} [Final time of simulation]
  step = 0.01 {caption: Step; category: Time; min: 0.01; max: 0.1; step: 0.001} [Time step of simlulation]

${CONTROL_EXPR.TOL}: 0.00005`;

/** Initial value problem templates */
export enum TEMPLATES {
  EMPTY = '',
  BASIC = TEMPLATE_BASIC,
  ADVANCED = TEMPLATE_ADVANCED,
  EXTENDED = TEMPLATE_EXTENDED,
};

/** Demo template */
export const DEMO_TEMPLATE = `${CONTROL_EXPR.NAME}: My model
${CONTROL_EXPR.DIF_EQ}:
  dx/dt = y * cos(p * t) + sin(t)
  dy/dt = x * sin(q * t) - exp(-t)

${CONTROL_EXPR.ARG}: t
  initial = 0 {min: 0; max: 3}
  final = 8 {min: 4; max: 20}
  step = 0.01 {min: 0.01; max: 0.1}

${CONTROL_EXPR.INITS}:
  x = 0 {min: 0; max: 5}
  y = 1 {min: 1; max: 7}

${CONTROL_EXPR.PARAMS}:
  p = -2 {min: -2; max: 2}
  q = -1 {min: -1; max: 1}`;
