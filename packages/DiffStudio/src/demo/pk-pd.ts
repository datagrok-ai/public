import {LINK} from '../ui-constants';
import {DemoModel} from './demo-model';

const MODEL = `#name: PK-PD
#equations:
  d(depot)/dt = -KA * depot
  d(centr)/dt = KA * depot - CL * C2 - Q * C2 + Q * C3
  d(peri)/dt  = Q * C2 - Q * C3
  d(eff)/dt  = Rate - Rate * (1 - C2/(EC50 + C2)) * eff

#expressions:
  C2 = centr / V2
  C3 = peri / V3

#output:
  t {caption: Time [h]}
  depot {caption: Depot}
  centr {caption: Central}
  peri {caption: Periferal}
  eff {caption: Effect}
  C2 {caption: Central concentration}
  C3 {caption: Peripheral concentration}

#loop:
  count = 10 {caption: count; category: Dosing; min: 1; max: 20} [Number of doses]
  depot += dose

#argument: t
  start = 0 {units: h; caption: begin; category: Misc; min: 0; max: 1} [Begin of dosing interval]
  final = 12 {units: h; caption: interval; category: Dosing; min: 5; max: 15} [End of dosing interval]
  step = 0.2 {units: h; caption: step; category: Misc; min: 0.1; max: 2} [Time step of simulation]  

#inits:  
  depot = 0 {category: Misc}
  centr = 0 {caption: central; category: Misc}
  peri = 0 {caption: peripheral; category: Misc}
  eff = 0.2 {caption: init effect; category: Misc}

#parameters:  
  dose = 1e4 {category: Dosing; min: 1e3; max: 2e4; step: 1e3} [Dosage]
  KA = 0.3 {caption: rate constant; category: PK parameters; min: 0.1; max: 1}
  CL = 2 {caption: clearance; category: PK parameters; min: 1; max: 5}
  V2 = 4 {caption: central volume; category: PK parameters; min: 1; max: 10} [Central compartment volume]
  Q = 1 {caption: inter rate; category: PK parameters; min: 0.1; max: 1} [Intercompartmental rate]
  V3 = 30 {caption: peri volume; category: PK parameters; min: 20; max: 40} [Peripheral compartment volume]
  EC50 = 8 {caption: effect; category: PD parameters; min: 1; max: 10}
  Rate = 0.2 {category: PD parameters; min: 0.1; max: 0.5} [Effective rate]`;

const UI_OPTS = {
  inputsTabDockRatio: 0.17,
  graphsDockRatio: 0.85,
};

const INFO = `# Model
Simulation of a two-compartment pharmacokinetic-pharmacodynamic
([PK-PD](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046)).
# Try
Interactive results based on input changes.
# Performance
Nonlinear systems of differential equations are solved within milliseconds.
# No-code
[Diff Studio](${LINK.DIF_STUDIO})
enables the creation of complex models without writing code.
# Learn more
* [Sensitivity analysis](${LINK.SENS_AN})
* [Parameter optimization](${LINK.FITTING})`;

export const PK_PD_DEMO = new DemoModel(MODEL, UI_OPTS, INFO);
