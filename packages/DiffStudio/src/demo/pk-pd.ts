/* eslint-disable max-len */
/** PK-PD demo model */

import {LINK} from '../ui-constants';
import {ModelInfo} from '../model';

/** PK-PD model specification */
const PK_PD_MODEL = `#name: PK-PD
#equations:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
  d(peri)/dt  = Q * C2 - Q * C3
  d(eff)/dt  = Rate - Rate * (1 - C2/(EC50 + C2)) * eff
  d(C2)/dt = (KA * depot - CL * C2 - Q * C2 + Q * C3) / V2
  d(C3)/dt = (Q * C2 - Q * C3) / V3

#output:
  t {caption: Time [h]}
  depot {caption: Depot}
  centr {caption: Central}
  peri {caption: Peripheral}
  eff {caption: Effect}
  C2 {caption: Central concentration}
  C3 {caption: Peripheral concentration}

#loop:
  count = 10 {caption: Count; category: Dosing; min: 1; max: 20} [Number of doses delivered in the loop]
  depot += dose

#argument: t
  start = 0  {units: h; caption: Begin;    category: Misc;   min: 0; max: 1}     [Start of the dosing interval]
  final = 12 {units: h; caption: Interval; category: Dosing; min: 5; max: 15}    [Length of one dosing interval (hours between doses)]
  step  = 0.2 {units: h; caption: Step;    category: Misc;   min: 0.1; max: 2}   [ODE solver time step]

#inits:  
  depot = 0   {caption: Depo; category: Misc}                     [Amount in the absorption site at t=0]
  centr = 0   {caption: Central;     category: Misc}              [Amount in the central compartment at t=0]
  peri  = 0   {caption: Peripheral;  category: Misc}              [Amount in the peripheral compartment at t=0]
  eff   = 0.2 {caption: Init effect; category: Misc}              [Initial effect / biomarker level]
  C2    = 0   {caption: Central;     category: Initial concentrations}  [Central concentration at t=0]
  C3    = 0   {caption: Peripheral;  category: Initial concentrations}  [Peripheral concentration at t=0]

#parameters:  
  dose = 1e4 {caption: Dose; category: Dosing; min: 1e3; max: 2e4; step: 1e3}             [Amount given at each dose]
  KA   = 0.3 {caption: Rate constant; category: PK parameters; min: 0.1; max: 1}          [Absorption rate constant (1/h)]
  CL   = 2   {caption: Clearance;     category: PK parameters; min: 1;   max: 5}          [Total clearance from central compartment]
  V2   = 4   {caption: Central volume; category: PK parameters; min: 1;   max: 10}        [Volume of the central compartment]
  Q    = 1   {caption: Inter rate;    category: PK parameters; min: 0.1; max: 1}          [Inter-compartmental clearance between central and peripheral]
  V3   = 30  {caption: Peri volume;   category: PK parameters; min: 20;  max: 40}         [Volume of the peripheral compartment]
  EC50 = 8   {caption: Effect;        category: PD parameters; min: 1;   max: 10}         [Drug concentration giving half-maximal effect]
  Rate = 0.2 {category: PD parameters; min: 0.1; max: 0.5}                                [Turnover rate of the effect variable]
  
#meta.solver: {method: 'lsoda'}`;

/** UI options */
const UI_OPTS = {
  inputsTabDockRatio: 0.17,
  graphsDockRatio: 0.85,
};

/** Info for help panel */
const INFO = `# Model
Simulation of a two-compartment pharmacokinetic-pharmacodynamic
([PK-PD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046)).
# Try
Interactive results based on input changes.
# Performance
Nonlinear systems of differential equations are solved within milliseconds.
# Quick reference

| Input | Meaning |
|---|---|
| Interval | Hours between doses |
| Dose | Amount per dose |
| Count | Number of doses |
| Central (concentration) | Central concentration at start |
| Peripheral (concentration) | Peripheral concentration at start |
| Rate constant | Drug absorption rate |
| Clearance | Drug elimination from blood |
| Central volume | Central compartment volume |
| Inter rate | Central-peripheral exchange rate |
| Peri volume | Peripheral compartment volume |
| Effect | Half-maximal effect concentration |
| Rate | Effect turnover rate |
| Begin | Dosing interval start |
| Step | ODE solver step size |
| Depo | Initial drug at absorption site |
| Central (amount) | Initial central compartment amount |
| Peripheral (amount) | Initial peripheral compartment amount |
| Init effect | Starting biomarker level |
`;

export const PK_PD_MODEL_INFO: ModelInfo = {
  equations: PK_PD_MODEL,
  uiOptions: UI_OPTS,
  info: INFO,
};
