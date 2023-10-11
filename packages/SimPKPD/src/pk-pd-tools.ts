/* PK-PD modelling tools.
 
  REFERENCES

    [1] Huixi Zou, et al, Application of Pharmacokinetic-Pharmacodynamic Modeling in Drug Delivery: 
        Development and Challenges (2020), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7348046/

    [2] RxODE, https://nlmixrdevelopment.github.io/RxODE/

    [3] W. Wang, et al, A Tutorial on RxODE: Simulating Differential Equation Pharmacometric Models in R (2016),
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728294/pdf/PSP4-5-03.pdf     
*/

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// Wasm-solving tools for PK-PD simulation
import { simulateOneCompartmentPK, simulateTwoCompartmentPK } from '../wasm/solving-tools';

const DEPOT_T = 'depot(t)';
const DEPOT = 'Depot';
const CENTR_T = 'centr(t)';
const CENTR = 'Central';
const PERI_T = 'peri(t)';
const PERI = 'Peripheral';
const EFF_T = 'eff(t)';
export const EFFECT = 'Effect';
const T_COL = 't';
export const TIME = 'Time [h]';
export const CENTR_CONC = 'Central concentration';
const PER_CONC = 'Peripheral concentration';
const POINTS_COUNTER = 100;

/** Compute concentrations */
function getConc(centrCol: DG.Column, periCol: DG.Column, V2: number, V3: number): {C2: DG.Column, C3: DG.Column} {
  if ((V2 === 0) || (V3 === 0))
    throw new Error('Zero volume: division by zero.');

  const centr = centrCol.getRawData();
  const peri = periCol.getRawData();
  const size = centr.length;

  if (size !== peri.length)
    throw new Error('The arrays peri & centr must be of the same size.');

  const C2 = new Float32Array(size);
  const C3 = new Float32Array(size);

  for (let i = 0; i < size; ++i) {
    C2[i] = centr[i] / V2;
    C3[i] = peri[i] / V3;
  }

  return {
    C2: DG.Column.fromFloat32Array(CENTR_CONC, C2),
    C3: DG.Column.fromFloat32Array(PER_CONC, C3)
  }
}

/** Get simulation function with respect to the model specified */
function getSimFn(compartments: string): any {
  switch (compartments) {
    case '1 compartment PK':
      return simulateOneCompartmentPK;

    case '2 compartment PK':
      return simulateTwoCompartmentPK;

    default:
      throw new Error('Unsupported PK/PD model.');
  }
}

/** Simulate PK-PD on the one dosage interval */
async function oneIntervalSimulation(compartments: string,
  initial: number, final: number, step: number,
  _depotInitial: number, _centrInitial: number, _periInitial: number, _effInitial: number, 
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, _KinVal: number, _KoutVal: number, _EC50Val: number): Promise<DG.DataFrame>
{
  // Specify modelling function
  const simFn = getSimFn(compartments);

  // Perform wasm-computations
  const res = await simFn(initial, final, step,
    _depotInitial, _centrInitial, _periInitial, _effInitial, 
    _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val) as DG.DataFrame;

  // Compute concentrations
  const conc = getConc(res.getCol(CENTR_T), res.getCol(PERI_T), _V2Val, _V3Val);

  // Add concentrations to the results
  res.columns.add(conc.C2);
  res.columns.add(conc.C3);

  return res;
}

/** Simulate PK-PD */
export async function simPKPD(compartments: string,
  dose: number, dosesCount: number, doseInterval: number,
  _KAVal: number, _CLVal: number, _V2Val: number, _QVal: number, _V3Val: number, _KinVal: number, _KoutVal: number, _EC50Val: number): Promise<DG.DataFrame>
{
  // Initial values
  let _depotInitial = dose;
  let _centrInitial = 0;
  let _periInitial = 0;
  let _effInitial = 1;

  // Independent variable range
  let initial = 0;
  let final = doseInterval;
  const step = doseInterval / POINTS_COUNTER;

  // Index of the last row
  let lastIdx;

  // Get simulation on the first interval
  const res = await oneIntervalSimulation(compartments, 
    initial, final, step,
    _depotInitial, _centrInitial, _periInitial, _effInitial, 
    _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val);

  // Get simulations on the next intervals
  for (let i = 1; i < dosesCount; ++i) {
    // Update independent variable range
    initial = final;
    final += doseInterval;

    lastIdx = res.rowCount - 1;
    
    // Update initial values
    _depotInitial = res.get(DEPOT_T, lastIdx) + dose;    
    _centrInitial = res.get(CENTR_T, lastIdx);
    _periInitial = res.get(PERI_T, lastIdx);
    _effInitial = res.get(EFF_T, lastIdx);

    // Compute simulation & update results
    res.append(
      await oneIntervalSimulation(compartments, 
        initial, final, step,
        _depotInitial, _centrInitial, _periInitial, _effInitial, 
        _KAVal, _CLVal, _V2Val, _QVal, _V3Val, _KinVal, _KoutVal, _EC50Val),
      true
    );
  }

  // Rename columns
  const timeCol = res.col(T_COL);
  timeCol!.name = TIME;
  
  const effCol = res.col(EFF_T);
  effCol!.name = EFFECT;

  const depotCol = res.col(DEPOT_T);
  depotCol!.name = DEPOT;

  const centrCol = res.col(CENTR_T);
  centrCol!.name = CENTR;

  const periCol = res.col(PERI_T);
  periCol!.name = PERI;

  return res;
}