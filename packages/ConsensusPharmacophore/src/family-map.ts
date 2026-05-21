/*
 * Single source of truth for the 7-family pharmacophore taxonomy.
 * Maps Datagrok 2D codes (from packages/Chem/files/pharmacophore-features.csv) to
 *   - rendered element symbol (Mol* applies CPK colors by element)
 *   - rendered resName (3-letter)
 *   - the 2D widget's hex color (kept here so the 3D legend matches the 2D widget)
 *   - canonical family name used in the consensus_model DataFrame
 *
 * Blueprint reference - Appendix A, Q16.
 */

export interface FamilyEntry {
  code: string;       // 'D' | 'A' | 'a' | 'H' | 'P' | 'N' | 'X'
  name: string;       // 'Donor' | 'Acceptor' | ...
  element: string;    // PDB element column (cols 77-78)
  resName: string;    // PDB resName (cols 18-20)
  hexColor: string;   // matches Chem widgets/pharmacophore-features.ts FAMILY_INFO
}

export const FAMILY_MAP: Record<string, FamilyEntry> = {
  D: {code: 'D', name: 'Donor',       element: 'N',  resName: 'HBD', hexColor: '#2196F3'},
  A: {code: 'A', name: 'Acceptor',    element: 'O',  resName: 'HBA', hexColor: '#E53935'},
  a: {code: 'a', name: 'Aromatic',    element: 'S',  resName: 'ARO', hexColor: '#4CAF50'},
  H: {code: 'H', name: 'Hydrophobic', element: 'C',  resName: 'HYD', hexColor: '#FFEB3B'},
  P: {code: 'P', name: 'Positive',    element: 'NA', resName: 'POS', hexColor: '#00BCD4'},
  N: {code: 'N', name: 'Negative',    element: 'CL', resName: 'NEG', hexColor: '#FF9800'},
  X: {code: 'X', name: 'Halogen',     element: 'BR', resName: 'HAL', hexColor: '#9C27B0'},
};

// Lookup by either the single-letter code OR the full name — Stage 5a emits names
// (consensus_model.family), Stage 4 internally uses codes.
export function resolveFamily(key: string): FamilyEntry {
  const trimmed = (key ?? '').toString().trim();
  if (FAMILY_MAP[trimmed]) return FAMILY_MAP[trimmed];
  for (const entry of Object.values(FAMILY_MAP))
    if (entry.name === trimmed) return entry;
  const low = trimmed.toLowerCase();
  for (const entry of Object.values(FAMILY_MAP))
    if (entry.name.toLowerCase() === low) return entry;
  if (low.includes('donor'))    return FAMILY_MAP.D;
  if (low.includes('accept'))   return FAMILY_MAP.A;
  if (low.includes('arom'))     return FAMILY_MAP.a;
  if (low.includes('hydroph'))  return FAMILY_MAP.H;
  if (low.includes('pos'))      return FAMILY_MAP.P;
  if (low.includes('neg'))      return FAMILY_MAP.N;
  if (low.includes('halo'))     return FAMILY_MAP.X;
  return FAMILY_MAP.H; // sensible default (Mol* draws C as grey, low visual noise)
}

export const FAMILY_CODES = ['D', 'A', 'a', 'H', 'P', 'N', 'X'] as const;
export const FAMILY_NAMES = FAMILY_CODES.map((c) => FAMILY_MAP[c].name);
