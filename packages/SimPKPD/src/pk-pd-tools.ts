/* PK-PD modelling tools. Further, we use notations from [3].
 
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

const DEPOT = 'Depot';
const CENTR = 'Central';
const PERI = 'Peripheral';
export const EFFECT = 'Effect';
export const TIME = 'Time [h]';
export const CENTR_CONC = 'Central concentration';
const PER_CONC = 'Peripheral concentration';
const POINTS_COUNTER = 100;

/** Simulation function type: notations from [3] are used.*/
type SimFn = (t0: number, t1: number, h: number, depot: number, centr: number, peri: number, eff: number, KA: number, CL: number, V2: number, Q: number, V3: number, Kin: number, Kout: number, EC50: number) => Promise<DG.DataFrame>;

/** Get concentrations that are evaluted outside wasm-computations.*/
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

/** Get simulation function with respect to the model specified. */
function getSimFn(compartments: string): SimFn {
  switch (compartments) {
    case '1 compartment PK':
      return simulateOneCompartmentPK;

    case '2 compartment PK':
      return simulateTwoCompartmentPK;

    default:
      throw new Error('Unsupported PK/PD model.');
  }
}

/** Simulate PK-PD on the one dosage interval. */
async function oneIntervalSimulation(compartments: string,
  t0: number, t1: number, h: number,
  depot: number, centr: number, peri: number, eff: number, 
  KA: number, CL: number, V2: number, Q: number, V3: number, Kin: number, Kout: number, E50: number): Promise<DG.DataFrame>
{
  // Specify modelling function
  const simFn = getSimFn(compartments);

  // Perform wasm-computations
  const res = await simFn(t0, t1, h, depot, centr, peri, eff, KA, CL, V2, Q, V3, Kin, Kout, E50) as DG.DataFrame;

  // Compute concentrations
  const conc = getConc(res.getCol(CENTR), res.getCol(PERI), V2, V3);

  // Add concentrations to the results
  res.columns.add(conc.C2);
  res.columns.add(conc.C3);

  return res;
}

/** Simulate PK-PD. */
export async function simPKPD(compartments: string,
  dose: number, dosesCount: number, doseInterval: number,
  KA: number, CL: number, V2: number, Q: number, V3: number, Kin: number, Kout: number, EC50: number): Promise<DG.DataFrame>
{
  // Initial values
  let depot = dose;
  let centr = 0;
  let peri = 0;
  let eff = 1;

  // Independent variable range
  let t0 = 0;
  let t1 = doseInterval;
  const h = doseInterval / POINTS_COUNTER;

  // Index of the last row
  let lastIdx = 0;

  // Get simulation on the first interval
  const res = await oneIntervalSimulation(compartments, t0, t1, h, depot, centr, peri, eff, KA, CL, V2, Q, V3, Kin, Kout, EC50);

  // Get simulations on the next intervals
  for (let i = 1; i < dosesCount; ++i) {
    // Update independent variable range
    t0 = t1;
    t1 += doseInterval;

    lastIdx = res.rowCount - 1;
    
    // Update initial values
    depot = res.get(DEPOT, lastIdx) + dose;    
    centr = res.get(CENTR, lastIdx);
    peri = res.get(PERI, lastIdx);
    eff = res.get(EFFECT, lastIdx);

    // Compute simulation & update results
    res.append(
      await oneIntervalSimulation(compartments, 
        t0, t1, h,
        depot, centr, peri, eff, 
        KA, CL, V2, Q, V3, Kin, Kout, EC50),
      true
    );
  }

  return res;
}