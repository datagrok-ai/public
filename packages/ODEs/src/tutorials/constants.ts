// Tutorials specific constants
import {CONTROL_EXPR} from '../constants';

/** Earth's population modeling */
export const POPULATION_MODEL = `${CONTROL_EXPR.NAME}: Population 
${CONTROL_EXPR.DIF_EQ}:
  dP/dt = k * P * (N - P)

${CONTROL_EXPR.ARG}: t
  initial = 2000 {category: Time; min: 0; max: 2010; units: year}
  final = 2100 {category: Time; min: 2010; max: 2200; units: year}
  step = 0.5 {category: Time; min: 0.5; max: 10; units: year}

${CONTROL_EXPR.INITS}:  
  P = 6.171 {caption: Population; category: Parameters; min: 5; max: 10; units: billion}

${CONTROL_EXPR.PARAMS}:
  k = 0.002 {caption: Growth rate; category: Parameters; min: 0.001; max: 0.01; units: 1/year}
  N = 12.5 {caption: Carrying capacity; category: Parameters; min: 2; max: 30; units: billion}`;
