/**
 * The Crux worker's molecule jobs (crux.worker.ts runs them; engine/molecule-jobs.ts and engine/worker-jobs.ts send
 * them): column-wide work that ran on the RDKit pool or on the main thread, in chunks across the Crux workers. Where
 * the RDKit worker has the code (conversions, CoordGen coordinates, flattening, beautifying) it runs on the Crux
 * module, which answers RDKit.js's calls; the rest is the main thread's own code (elemental analysis, the R-group
 * dialog's Core column, the SAR matrix's atom counts) moved here as it stands. Each job gives its results and the rows
 * the RDKit path does (a molecule Crux does not read, a call it refuses).
 */
import type {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getMolSafe} from '../../utils/mol-creation_rdkit';
import {hasNewLines} from '../../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../../utils/chem-constants';
import {V2000_ATOM_NAME_LEN, V2000_ATOM_NAME_POS} from '../../constants';
import {isV3000Block, v3000Atoms} from '../../utils/molblock-atoms';
import {notConverted} from '../../rdkit-service/rdkit-service-worker-substructure';
import {WorkerJob, throwIfTrap} from '../crux-api';

/** A job's results, one per molecule, and the rows the RDKit path does. */
export interface MoleculeJobResult<T> {results: T[], rest: number[]}

/** Each molecule's element counts and charge (utils/elemental-analysis-utils.ts `getAtomsColumn`'s), for a chunk. */
export interface ElementCounts {
  /** Per element symbol (as the V2000 MOL block writes it, `R` included), its atoms per molecule. */
  elements: {[symbol: string]: Int32Array};
  /** The sum of the atoms' formal charges per molecule. */
  charges: Int32Array;
}

/** The rows of `results` with a molecule and no result: the RDKit pool's code gives '' where it reads no molecule or
 * its call throws, and the RDKit path does them. */
function emptyRows(molecules: string[], results: string[]): number[] {
  return molecules.flatMap((m, i) => m && !results[i] ? [i] : []);
}

/**
 * The atoms of a V2000 MOL block by symbol and their charge, as `getAtomsColumn` reads RDKit's MOL block: the counts
 * line's atom count, each atom line's symbol; the charge OpenChemLib reads from it (`M  CHG`, where RDKit writes every
 * charge; the atom lines' charge field, 0 in RDKit's blocks, otherwise).
 */
function countAtoms(block: string, row: number, out: ElementCounts, length: number): void {
  if (isV3000Block(block)) { // a molecule V2000 cannot hold (enhanced stereo, > 999 atoms), as getAtomsColumn counts it
    const {symbols, charge} = v3000Atoms(block);
    for (const symbol of symbols) {
      out.elements[symbol] ??= new Int32Array(length);
      ++out.elements[symbol][row];
    }
    out.charges[row] = charge;
    return;
  }
  let curPos = 0;
  curPos = block.indexOf('\n', curPos) + 1;
  curPos = block.indexOf('\n', curPos) + 1;
  curPos = block.indexOf('\n', curPos) + 1;
  const atomCounts = parseInt(block.substring(curPos, curPos + 3));
  const atomCharges: number[] = [];
  for (let atomRowI = 0; atomRowI < atomCounts; atomRowI++) {
    curPos = block.indexOf('\n', curPos) + 1;
    const symbol = block.substring(curPos + V2000_ATOM_NAME_POS, curPos + V2000_ATOM_NAME_POS + V2000_ATOM_NAME_LEN)
      .trim();
    out.elements[symbol] ??= new Int32Array(length);
    ++out.elements[symbol][row];
    const field = parseInt(block.substring(curPos + 36, curPos + 39));
    atomCharges.push(field > 0 && field < 8 && field !== 4 ? 4 - field : 0);
  }
  const chg = [...block.matchAll(/^M {2}CHG {1,2}\d+((?: +-?\d+ +-?\d+)+)/gm)];
  if (chg.length) {
    atomCharges.fill(0);
    for (const line of chg) {
      const values = line[1].trim().split(/\s+/).map((v) => parseInt(v));
      for (let k = 0; k + 1 < values.length; k += 2) {
        if (values[k] >= 1 && values[k] <= atomCounts)
          atomCharges[values[k] - 1] = values[k + 1];
      }
    }
  }
  out.charges[row] = atomCharges.reduce((a, b) => a + b, 0);
}

/** The jobs by name (the Crux worker's `op`). */
export const MOLECULE_JOBS: {[op: string]: WorkerJob} = {
  /**
   * The column conversions (engine/worker-jobs.ts `convertMolNotationByEngine`: RdKitServiceWorker.convertMolNotation),
   * the rows left unconverted (a molecule Crux cannot read, a writer it refuses) and what such a row holds.
   */
  convertMolNotation: async ({rdkitWorker}, molecules: string[], targetNotation: string, kekulize: boolean):
    Promise<{results: string[], none: string, unconverted: number[]}> => {
    const results = await rdkitWorker().convertMolNotation(molecules, targetNotation, kekulize);
    const none = notConverted(targetNotation);
    return {results, none, unconverted: molecules.flatMap((m, i) => m && results[i] === none ? [i] : [])};
  },
  /** Transform | Recalculate Coordinates by CoordGen (RdKitServiceWorker.getCoordGenCoords: `get_new_coords(true)`). */
  coordGenCoords: async ({rdkitWorker}, molecules: string[]): Promise<MoleculeJobResult<string>> => {
    const results = await rdkitWorker().getCoordGenCoords(molecules);
    return {results, rest: emptyRows(molecules, results)};
  },
  /** Transform | Flatten Molecules (RdKitServiceWorker.flattenMolecules: SMILES without stereo). */
  flattenMolecules: async ({rdkitWorker}, molecules: string[]): Promise<MoleculeJobResult<string>> => {
    const results = await rdkitWorker().flattenMolecules(molecules);
    return {results, rest: emptyRows(molecules, results)};
  },
  /** `beautifyMols` (RdKitServiceWorker.beautifyMoleculesV3K: new coordinates, normalized and straightened, V3000). */
  beautifyMols: async ({rdkitWorker}, molecules: string[]): Promise<MoleculeJobResult<string>> => {
    const results = await rdkitWorker().beautifyMoleculesV3K(molecules);
    return {results, rest: emptyRows(molecules, results)};
  },
  /**
   * Elemental analysis (`getAtomsColumn`): each molecule read as getMolSafe reads it, its Hs added, its V2000 MOL block
   * written and counted ({@link countAtoms}); `rest`, the molecules Crux does not read or write. An empty one counts
   * nothing (RDKit's empty molecule).
   */
  elementalAnalysis: ({module}, molecules: string[]): ElementCounts & {rest: number[]} => {
    const out: ElementCounts = {elements: {}, charges: new Int32Array(molecules.length)};
    const rest: number[] = [];
    molecules.forEach((molecule, row) => {
      if (!molecule)
        return;
      const mol = getMolSafe(molecule, {}, module).mol;
      if (!mol) {
        rest.push(row);
        return;
      }
      try {
        try {
          mol.add_hs_in_place();
        } catch (e) {
          throwIfTrap(e);
        }
        countAtoms(mol.get_molblock(), row, out, molecules.length);
      } catch (e) {
        throwIfTrap(e);
        rest.push(row);
      } finally {
        mol.delete();
      }
    });
    return {...out, rest};
  },
  /**
   * The R-group dialog's Core column (analysis/r-group-analysis.ts): each core's MOL block (`get_mol`, else `get_qmol`)
   * with its ISO line written as RGP; '' for an empty core, null for a SMILES over Chem's length limit (no molecule);
   * `rest`, the cores Crux reads no molecule from.
   */
  coreMolblocks: ({module}, cores: string[]): MoleculeJobResult<string | null> => {
    const rest: number[] = [];
    const results = cores.map((core, row) => {
      if (!core)
        return '';
      if (!hasNewLines(core) && core.length > MAX_SMILES_LENGTH)
        return null;
      let mol: RDMol | null = null;
      try {
        mol = module.get_mol(core);
        if (!mol)
          mol = module.get_qmol(core);
        if (mol)
          return mol.get_molblock().replace('ISO', 'RGP');
      } catch (e) {
        throwIfTrap(e);
      } finally {
        mol?.delete();
      }
      rest.push(row);
      return '';
    });
    return {results, rest};
  },
  /**
   * The SAR matrix's heavy-atom counts (analysis/sar-matrix `cachedAtomCount`): a valid molecule's atoms; `rest`, the
   * structures Crux reads none from (counted 0 by RDKit too where it reads none).
   */
  atomCounts: ({module}, structures: string[]): MoleculeJobResult<number> => {
    const rest: number[] = [];
    const results = structures.map((s, row) => {
      if (s && !hasNewLines(s) && s.length > MAX_SMILES_LENGTH) { // (the wasm could run out of memory): RDKit's
        rest.push(row);
        return 0;
      }
      let mol: RDMol | null = null;
      try {
        mol = module.get_mol(s);
        const count = mol && mol.is_valid() ? mol.get_num_atoms() : 0;
        if (count === 0 && s)
          rest.push(row);
        return count;
      } catch (e) {
        throwIfTrap(e);
        rest.push(row);
        return 0;
      } finally {
        mol?.delete();
      }
    });
    return {results, rest};
  },
};
