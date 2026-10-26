/**
 * The Crux worker's calls of chem-wasm's own calculations (crux.worker.ts runs them, CruxService sends them): the
 * descriptors as the chem-chem container gives them, the scaffold tree, Chem's property table, QED, FilterCatalog
 * alerts, InChI's input and 3D coordinates.
 */
import type {JSMol} from '../../../vendor/chem-wasm/chem_wasm.js';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe} from '../../utils/mol-creation_rdkit';
import {DescriptorValues, INT_NONE, PackedDescriptorValues, WorkerJob, throwIfTrap} from '../crux-api';

/** The container-shaped descriptor values as typed arrays (the main thread gets their buffers, not copies). */
function packDescriptors(values: {[name: string]: DescriptorValues}, transfer: Transferable[]):
  {[name: string]: PackedDescriptorValues} {
  const packed: {[name: string]: PackedDescriptorValues} = {};
  for (const [name, {type, value}] of Object.entries(values)) {
    const array = type === 'int' ? new Int32Array(value.length) : new Float64Array(value.length);
    const none = type === 'int' ? INT_NONE : NaN;
    for (let i = 0; i < value.length; i++)
      array[i] = value[i] ?? none;
    packed[name] = {type, value: array};
    transfer.push(array.buffer);
  }
  return packed;
}

/**
 * InChI's input of each molecule, as RDKit's MolToInchi fills it (chem-wasm's `get_inchi_input`) for crux-sketch's
 * InChI module to run the library on (engine/worker-jobs.ts); '' where RDKit's `get_inchi` gives '' (an empty cell, no
 * input, a throw) and null where Crux cannot read the molecule. The inputs' buffers are transferred.
 */
function inchiInputs(module: RDModule, molecules: string[], transfer: Transferable[]): (Uint8Array | '' | null)[] {
  return molecules.map((m) => {
    const mol = m ? getMolSafe(m, {}, module).mol : null;
    if (!mol)
      return m ? null : '';
    try {
      const input = (mol as unknown as JSMol).get_inchi_input('');
      if (input)
        transfer.push(input.buffer);
      return input ?? '';
    } catch (e) {
      throwIfTrap(e);
      return '';
    } finally {
      mol.delete();
    }
  });
}

/** The jobs by name (the Crux worker's `op`). */
export const CALCULATION_JOBS: {[op: string]: WorkerJob} = {
  /** The container's descriptors (descriptors/descriptors-engine.ts), packed. */
  'calcDescriptors': ({crux, transfer}, molecules: (string | null)[], names: string[]) =>
    packDescriptors(JSON.parse(crux.calcDescriptors(molecules, names)), transfer),
  /** Diagnosis: {@link CALCULATION_JOBS.calcDescriptors} with the Crux call and the parse (and packing) timed apart. */
  'calcDescriptors::timed': ({crux, transfer}, molecules: (string | null)[], names: string[]) => {
    const started = performance.now();
    const text = crux.calcDescriptors(molecules, names);
    const computed = performance.now();
    return {values: packDescriptors(JSON.parse(text), transfer), wasmMs: computed - started,
      parseMs: performance.now() - computed};
  },
  /** scripts/scaffold_tree.py's JSON text. */
  'scaffoldTree': ({crux}, smiles: string[], ringCutoff: number, discharge: boolean): string =>
    crux.scaffoldTree(smiles, ringCutoff, discharge),
  /** Chem's property table (widgets/properties.ts). */
  'chemProperties': ({crux}, molecules: (string | null)[], names: string[]) =>
    JSON.parse(crux.chemProperties(molecules, names)),
  /** Drug-likeness: QED and its properties (widgets/drug-likeness.ts). */
  'chemQed': ({crux}, molecules: (string | null)[]) => JSON.parse(crux.chemQed(molecules)),
  /** Toxicity: FilterCatalog's Brenk and NIH alerts (widgets/toxicity.ts). */
  'chemAlerts': ({crux}, molecules: (string | null)[], catalogs: string[], atoms: boolean) =>
    JSON.parse(crux.chemAlerts(molecules, catalogs, atoms)),
  /** InChI's input per molecule (engine/worker-jobs.ts `inchisByEngine`). */
  'inchiInputs': ({module, transfer}, molecules: string[]) => inchiInputs(module, molecules, transfer),
  /** The 3D Structure pane's coordinates (widgets/structure3d.ts); a timeout throws. */
  'chemEmbed3d': ({crux}, molecule: string, options: string): string | null =>
    crux.chemEmbed3d(molecule, options) ?? null,
};
