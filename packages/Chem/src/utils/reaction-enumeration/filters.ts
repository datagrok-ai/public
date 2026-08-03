import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ProductsSpecs} from './config';

const NONMETAL_Z = new Set([1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54, 85, 86]);
const HALOGEN_Z = new Set([9, 17, 35, 53, 85]);

const Z_TO_SYMBOL: Record<number, string> = {
  1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
  11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
  19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
  29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
  37: 'Rb', 38: 'Sr', 47: 'Ag', 48: 'Cd', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe',
  55: 'Cs', 56: 'Ba', 78: 'Pt', 79: 'Au', 80: 'Hg', 82: 'Pb', 83: 'Bi', 85: 'At', 86: 'Rn',
};

function isHetero(z: number): boolean { return z !== 1 && z !== 6; }
function isMetal(z: number): boolean { return !NONMETAL_Z.has(z); }
function isHalogen(z: number): boolean { return HALOGEN_Z.has(z); }

export interface MolStats {
  numHeavy: number;
  countByZ: Map<number, number>;
  numAromatic: number;
  numUnsaturatedNonAromatic: number;
  hasIsotope: boolean;
  hasCharge: boolean;
  hasRadical: boolean;
  symbols: string[];
}

export function computeMolStats(mol: RDMol): MolStats {
  let raw: string;
  try { raw = mol.get_json(); } catch (e) {
    throw new Error(`get_json failed: ${e instanceof Error ? e.message : String(e)}`);
  }
  let json: any;
  try { json = JSON.parse(raw); } catch (e) {
    throw new Error(`get_json returned invalid JSON: ${e instanceof Error ? e.message : String(e)}`);
  }
  const m = json.molecules?.[0];
  if (!m) throw new Error('get_json returned no molecule');

  const defaults = json.defaults?.atom ?? {};
  const dz = defaults.z ?? 6;
  const dchg = defaults.chg ?? 0;
  const diso = defaults.isotope ?? 0;
  const drad = defaults.nRad ?? 0;
  const dbo = json.defaults?.bond?.bo ?? 1;

  const aromaticAtoms = new Set<number>();
  const aromaticBonds = new Set<number>();
  for (const e of (m.extensions ?? []) as any[]) {
    if (e.name === 'rdkitRepresentation') {
      for (const i of e.aromaticAtoms ?? []) aromaticAtoms.add(i);
      for (const i of e.aromaticBonds ?? []) aromaticBonds.add(i);
    }
  }

  const stats: MolStats = {
    numHeavy: 0, countByZ: new Map(),
    numAromatic: 0, numUnsaturatedNonAromatic: 0,
    hasIsotope: false, hasCharge: false, hasRadical: false,
    symbols: [],
  };

  const atoms = (m.atoms ?? []) as any[];
  for (let i = 0; i < atoms.length; i++) {
    const a = atoms[i];
    const z = a.z ?? dz;
    const chg = a.chg ?? dchg;
    const iso = a.isotope ?? diso;
    const rad = a.nRad ?? drad;
    if (z !== 1) stats.numHeavy++;
    stats.countByZ.set(z, (stats.countByZ.get(z) ?? 0) + 1);
    stats.symbols.push(Z_TO_SYMBOL[z] ?? `?${z}`);
    if (iso !== 0) stats.hasIsotope = true;
    if (chg !== 0) stats.hasCharge = true;
    if (rad !== 0) stats.hasRadical = true;
    if (aromaticAtoms.has(i)) stats.numAromatic++;
  }

  const bonds = (m.bonds ?? []) as any[];
  for (let i = 0; i < bonds.length; i++) {
    const bo = bonds[i].bo ?? dbo;
    if (bo > 1 && !aromaticBonds.has(i)) stats.numUnsaturatedNonAromatic++;
  }

  return stats;
}

export interface FilterResult { pass: boolean; reason?: string; }

export function applyProductFilters(
  stats: MolStats,
  specs: ProductsSpecs,
  exclusionQmols: RDMol[],
  mol: RDMol,
): FilterResult {
  if (specs.max_num_heavy_atoms >= 0 && stats.numHeavy > specs.max_num_heavy_atoms)
    return {pass: false, reason: 'max_num_heavy_atoms'};

  const carbonCount = stats.countByZ.get(6) ?? 0;
  if (specs.min_num_carbon_atoms >= 0 && carbonCount < specs.min_num_carbon_atoms)
    return {pass: false, reason: 'min_num_carbon_atoms'};
  if (specs.max_num_carbon_atoms >= 0 && carbonCount > specs.max_num_carbon_atoms)
    return {pass: false, reason: 'max_num_carbon_atoms'};

  let hetero = 0; let halogen = 0; let metal = 0;
  for (const [z, n] of stats.countByZ) {
    if (isHetero(z)) hetero += n;
    if (isHalogen(z)) halogen += n;
    if (isMetal(z)) metal += n;
  }
  if (specs.max_num_hetero_atoms >= 0 && hetero > specs.max_num_hetero_atoms)
    return {pass: false, reason: 'max_num_hetero_atoms'};
  if (specs.max_num_nitrogen >= 0 && (stats.countByZ.get(7) ?? 0) > specs.max_num_nitrogen)
    return {pass: false, reason: 'max_num_nitrogen'};
  if (specs.max_num_sulfur >= 0 && (stats.countByZ.get(16) ?? 0) > specs.max_num_sulfur)
    return {pass: false, reason: 'max_num_sulfur'};
  if (specs.max_num_oxygen >= 0 && (stats.countByZ.get(8) ?? 0) > specs.max_num_oxygen)
    return {pass: false, reason: 'max_num_oxygen'};
  if (specs.max_num_metals >= 0 && metal > specs.max_num_metals)
    return {pass: false, reason: 'max_num_metals'};
  if (specs.max_num_halogens >= 0 && halogen > specs.max_num_halogens)
    return {pass: false, reason: 'max_num_halogens'};
  if (specs.max_num_aromatic_atoms >= 0 && stats.numAromatic > specs.max_num_aromatic_atoms)
    return {pass: false, reason: 'max_num_aromatic_atoms'};
  if (specs.max_num_unsaturated_nonaromatic_bonds >= 0 &&
      stats.numUnsaturatedNonAromatic > specs.max_num_unsaturated_nonaromatic_bonds)
    return {pass: false, reason: 'max_num_unsaturated_nonaromatic_bonds'};

  if (specs.only_these_atoms_allowed && specs.only_these_atoms_allowed.length > 0) {
    const allowed = new Set(specs.only_these_atoms_allowed.map((s) => s.trim()));
    for (const sym of stats.symbols) {
      if (!allowed.has(sym)) return {pass: false, reason: 'only_these_atoms_allowed'};
    }
  }

  if (specs.remove_radicals && stats.hasRadical)
    return {pass: false, reason: 'remove_radicals'};
  if (specs.remove_charged_species && stats.hasCharge)
    return {pass: false, reason: 'remove_charged_species'};

  for (const qmol of exclusionQmols) {
    let match: string;
    try { match = mol.get_substruct_match(qmol); } catch { continue; }
    if (match && match !== '' && match !== '{}') {
      try {
        const obj = JSON.parse(match);
        if (obj && Array.isArray(obj.atoms) && obj.atoms.length > 0)
          return {pass: false, reason: 'exclusion_smarts'};
      } catch {
        if (match.length > 2) return {pass: false, reason: 'exclusion_smarts'};
      }
    }
  }

  return {pass: true};
}

export function stripIsotopesFromSmiles(smiles: string): string {
  return smiles.replace(/\[(\d+)([A-Za-z])/g, '[$2');
}

export function buildExclusionQmols(rdkit: RDModule, smartsList: string[]): {qmols: RDMol[]; dispose: () => void} {
  const qmols: RDMol[] = [];
  const sources: string[] = [];
  for (const s of smartsList) {
    const trimmed = (s ?? '').trim();
    if (!trimmed) continue;
    try {
      const q = rdkit.get_qmol(trimmed);
      if (q && q.is_valid()) { qmols.push(q); sources.push(trimmed); }
      else q?.delete();
    } catch {
      // skip invalid SMARTS
    }
  }
  let disposed = false;
  return {
    qmols,
    dispose: () => {
      if (disposed) return;
      disposed = true;
      for (let i = 0; i < qmols.length; i++) {
        try { qmols[i].delete(); }
        catch (e) {
          console.warn(`exclusion qmol[${i}] delete failed (SMARTS="${sources[i]}"): ` +
            `${e instanceof Error ? e.message : String(e)}`);
        }
      }
      qmols.length = 0;
      sources.length = 0;
    },
  };
}
