/* eslint-disable max-len */
import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {hasNewLines, isMolBlock} from './chem-common';
import {MAX_SMILES_LENGTH} from './chem-constants';
import {MolfileHandler} from '@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler';

export interface IMolContext {
  mol: RDMol | null; // null when molString is invalid
  kekulize: boolean;
  isQMol: boolean;
  useMolBlockWedging: boolean;
}

export function isFragment(molString: string) {
  if (isMolBlock(molString))
    return MolfileHandler.getInstance(molString).isFragment();
  else
    return !!molString.match(/\[.?:|\*.?\]/g);
}

/** A CXSMILES extension (` |...|` after the SMILES): its atom labels hold `$` ($Q_e$) and its stereo groups `&` (&1:2),
 * which are no SMARTS. */
const CX_EXTENSION = /\s+\|[^|]*\|\s*$/;

/** A molecule string without its CXSMILES extension, if it has one. */
export function withoutCxExtension(molString: string): string {
  return molString.replace(CX_EXTENSION, '');
}

/** Whether a molecule string is a query: a query molblock, or a string with SMARTS syntax (a `[#n]` atom, a recursive
 * `$()`, the logical `&`, `;`, `,` and `!`, the any bond `~`) outside a CXSMILES extension. */
export function _isSmarts(molString: string): boolean {
  if (isMolBlock(molString))
    return MolfileHandler.getInstance(molString).isQuery();
  else
    return !!withoutCxExtension(molString).match(/\[.?#\d|\$|&|;|,|!.?]|~/g);
}

/** The element symbols RDKit's SMILES reader takes in a bracket atom (H to Og), and the aromatic ones. */
const SMILES_ELEMENT = 'H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|' +
  'Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|' +
  'W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|' +
  'Cn|Nh|Fl|Mc|Lv|Ts|Og|as|se|si|te|[bcnops*]';

/** A string of SMILES tokens alone: atoms of the organic subset and `*`; bracket atoms in SMILES' order (isotope, an
 * element, `@` or `@@`, H count, charge, class); a bond or a dative bond not followed by another (two in a row, `=\` or
 * `-/`, are SMARTS' conjunction); ring closures, branches and dots. */
const SMILES_SYNTAX = new RegExp(`^(?:\\[\\d{0,3}(?:${SMILES_ELEMENT})@{0,2}(?:H\\d{0,2})?(?:[+-]\\d{0,2}|\\+\\+|--)?(?::\\d+)?\\]` +
  '|Cl|Br|[BCNOPSFIbcnops*]|(?:->|<-|[-=#:/\\\\])(?![-=#:/\\\\<])|[.()]|%\\(\\d+\\)|%\\d\\d|\\d)*$');

/** Whether a string is written in SMILES syntax alone: no SMARTS primitive or operator, no CXSMILES extension, no name.
 * RDKit reads such a string as a molecule or not at all, never as a query only; a SMARTS that no syntax tells from SMILES
 * holds a token outside it ([NX3], [R], A, a, the ring bond @, an H count after a charge). */
export function isSmilesSyntax(s: string): boolean {
  return SMILES_SYNTAX.test(s);
}

/** Whether RDKit reads `s` as a query and not as a molecule, even unsanitized: a SMARTS of primitives alone, which no
 * syntax tells from SMILES (Daylight's [NX3][CX3]=[OX1], [R], A-A, C@C). */
export function readsAsQueryOnly(s: string, rdKitModule: RDModule): boolean {
  if (!s || isMolBlock(s) || (!hasNewLines(s) && s.length > MAX_SMILES_LENGTH))
    return false;
  let mol: RDMol | null = null;
  try {
    mol = rdKitModule.get_mol(s, JSON.stringify({sanitize: false}));
  } catch {
    mol = null;
  }
  if (mol !== null) {
    mol.delete();
    return false;
  }
  let query: RDMol | null = null;
  try {
    query = rdKitModule.get_qmol(s);
  } catch {
    query = null;
  }
  query?.delete();
  return query !== null;
}

export function hasRadicals(mol: RDMol): boolean {
  try {
    return /"nRad":[1-9]/.test(mol.get_json());
  } catch {
    return false;
  }
}

export function getMolSafe(molString: string, details: object = {}, rdKitModule: RDModule,
  warnOff: boolean = true, onError?: (e: unknown) => void): IMolContext {
  // `onError` sees what the module threw (Chem's Crux module: a trap marks Crux failed); the molecule is null then
  if (molString && !hasNewLines(molString) && molString.length > MAX_SMILES_LENGTH)
    return {mol: null, kekulize: true, isQMol: false, useMolBlockWedging: false}; // do not attempt to parse very long SMILES, will cause MOB.
  let isQMol = false;
  const kekulizeProp = (details as any).kekulize;
  let kekulize: boolean = typeof kekulizeProp === 'boolean' ? kekulizeProp : true;
  let useMolBlockWedging: boolean = false;
  let mol: RDMol | null = null;

  try {
    mol = rdKitModule.get_mol(molString, JSON.stringify(details));
  } catch (e) {
    onError?.(e);
  }
  if (!mol && kekulize) {
    kekulize = false; //Pyrrole cycles
    try {
      mol = rdKitModule.get_mol(molString, JSON.stringify({...details, kekulize}));
    } catch (e) {
      onError?.(e);
    }
  }
  if (!mol) {
    kekulize = true; //previous step failed, so retrun kekulize parameter to true
    isQMol = true;
    try {
      mol = rdKitModule.get_qmol(molString);
      isQMol = true;
    } catch (e) {
      onError?.(e);
    }
  }
  if (mol)
    useMolBlockWedging = (mol.has_coords() === 2);
  else if (!warnOff)
    console.error('Chem | In getMolSafe: RDKit.get_mol crashes on a molString: `' + molString + '`');
  return {mol, kekulize, isQMol, useMolBlockWedging};
}


/** Whether an atom of a molblock carries an H count query, as RDKit's readers take one: a V2000 atom line's hhh field
 * (columns 43-45), a V3000 atom's HCOUNT (0 or blank: none). */
function hasHCountQuery(molBlock: string): boolean {
  const lines = molBlock.split(/\r?\n/);
  const counts = lines[3] ?? '';
  if (counts.includes('V3000')) {
    let inAtoms = false;
    for (let i = 4; i < lines.length; i++) {
      let line = lines[i];
      if (!line.startsWith('M  V30 '))
        continue;
      while (line.endsWith('-') && i + 1 < lines.length && lines[i + 1].startsWith('M  V30 '))
        line = line.slice(0, -1) + lines[++i].slice(7);
      const body = line.slice(7).trim();
      if (body === 'BEGIN ATOM')
        inAtoms = true;
      else if (body === 'END ATOM')
        return false;
      else if (inAtoms) {
        const hCount = body.match(/(?:^|\s)HCOUNT=(-?\d+)/);
        if (hCount && parseInt(hCount[1]) !== 0)
          return true;
      }
    }
    return false;
  }
  const atomCount = parseInt(counts.substring(0, 3));
  for (let i = 4; i < 4 + (atomCount || 0) && i < lines.length; i++) {
    const hhh = lines[i].length >= 45 ? lines[i].substring(42, 45).trim() : '';
    if (hhh !== '' && parseInt(hhh) !== 0)
      return true;
  }
  return false;
}

/** V2000's atom block charges (ccc), by their code. */
const V2000_CHARGES: {[ccc: number]: number} = {1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3};
/** A bond's order by its molblock type (an aromatic bond 1.5); a query or other bond counts 1. */
const BOND_ORDERS: {[type: number]: number} = {1: 1, 2: 2, 3: 3, 4: 1.5};

interface MolBlockGraph {
  /** by their number in the molblock: the symbol, the charge, the atom's line (V2000) */
  atoms: Map<number, {symbol: string, charge: number, line: number}>;
  /** [atom, atom, type] */
  bonds: [number, number, number][];
}

/** A molblock's atoms and bonds, read from its lines as far as the checks below need them (V2000: an atom's charge is its
 * ccc unless M  CHG lines give the charges; V3000: CHG=). */
function molBlockGraph(lines: string[]): MolBlockGraph {
  const atoms: MolBlockGraph['atoms'] = new Map();
  const bonds: MolBlockGraph['bonds'] = [];
  const counts = lines[3] ?? '';
  if (counts.includes('V3000')) {
    let block = '';
    for (let i = 4; i < lines.length; i++) {
      let line = lines[i];
      if (!line.startsWith('M  V30 '))
        continue;
      const first = i;
      while (line.endsWith('-') && i + 1 < lines.length && lines[i + 1].startsWith('M  V30 '))
        line = line.slice(0, -1) + lines[++i].slice(7);
      const body = line.slice(7).trim();
      if (body.startsWith('BEGIN ') || body.startsWith('END ')) {
        block = body.startsWith('BEGIN ') ? body.slice(6) : '';
        continue;
      }
      const fields = body.split(/\s+/);
      if (block === 'ATOM') {
        const charge = body.match(/(?:^|\s)CHG=(-?\d+)/);
        atoms.set(parseInt(fields[0]), {symbol: fields[1], charge: charge ? parseInt(charge[1]) : 0, line: first});
      } else if (block === 'BOND')
        bonds.push([parseInt(fields[2]), parseInt(fields[3]), parseInt(fields[1])]);
    }
    return {atoms, bonds};
  }
  const atomCount = parseInt(counts.substring(0, 3)) || 0;
  const bondCount = parseInt(counts.substring(3, 6)) || 0;
  for (let n = 1; n <= atomCount && 3 + n < lines.length; n++) {
    const line = lines[3 + n];
    atoms.set(n, {symbol: line.substring(31, 34).trim(), charge: V2000_CHARGES[parseInt(line.substring(36, 39))] ?? 0,
      line: 3 + n});
  }
  for (let i = 4 + atomCount; i < 4 + atomCount + bondCount && i < lines.length; i++)
    bonds.push([parseInt(lines[i].substring(0, 3)), parseInt(lines[i].substring(3, 6)), parseInt(lines[i].substring(6, 9))]);
  let chargeLines = false;
  for (let i = 4 + atomCount + bondCount; i < lines.length; i++) {
    if (!lines[i].startsWith('M  CHG'))
      continue;
    if (!chargeLines) { // the first M  CHG line resets the atom block's charges
      atoms.forEach((atom) => atom.charge = 0);
      chargeLines = true;
    }
    const fields = lines[i].substring(6).trim().split(/\s+/).map((field) => parseInt(field));
    for (let k = 0; k < fields[0] && 2 + 2 * k < fields.length; k++) {
      const atom = atoms.get(fields[1 + 2 * k]);
      if (atom)
        atom.charge = fields[2 + 2 * k];
    }
  }
  return {atoms, bonds};
}

/** Whether a molblock (V2000 or V3000) holds an uncharged N with five bonds' worth of valence: a nitro drawn N(=O)=O, an
 * azide drawn N=N#N, a pyridine N-oxide drawn N=O. */
function hasPentavalentN(molBlock: string): boolean {
  const {atoms, bonds} = molBlockGraph(molBlock.split(/\r?\n/));
  const valence = new Map<number, number>();
  for (const [a, b, type] of bonds) {
    valence.set(a, (valence.get(a) ?? 0) + (BOND_ORDERS[type] ?? 1));
    valence.set(b, (valence.get(b) ?? 0) + (BOND_ORDERS[type] ?? 1));
  }
  for (const [n, atom] of atoms) {
    if (atom.symbol === 'N' && atom.charge === 0 && (valence.get(n) ?? 0) >= 5)
      return true;
  }
  return false;
}

/** A V2000 molblock with each uncharged atom that carries an H count query (hhh: at least hhh - 1 H) given that many H by
 * its valence field (vvv), where that field is blank and the atom's bonds are all single, double or triple. RDKit's V2000
 * reader gives an atom with an hhh no implicit H, so that sanitizing made a ring NH drawn Kekulé (OpenChemLib's indole) a
 * radical whose ring was not aromatic: the query matched nothing (its V3000 twin, HCOUNT, matched the indoles). The H
 * count query itself is RDKit's as before. A charged atom is left: RDKit's query of an atom with an hhh holds no charge,
 * so that a pyridinium NH+ read aromatic matched the pyridones' neutral NH. Any other molblock is returned as it is. */
function withHCountHs(molBlock: string): string {
  const lines = molBlock.split(/\r?\n/);
  if ((lines[3] ?? '').includes('V3000'))
    return molBlock;
  const {atoms, bonds} = molBlockGraph(lines);
  let changed = false;
  for (const [n, atom] of atoms) {
    const line = lines[atom.line];
    const hhh = line.length >= 45 ? parseInt(line.substring(42, 45)) : NaN;
    if (!(hhh >= 2) || (parseInt(line.substring(48, 51)) || 0) !== 0 || atom.charge !== 0)
      continue;
    const own = bonds.filter(([a, b]) => a === n || b === n);
    if (own.some(([, , type]) => !(type >= 1 && type <= 3)))
      continue;
    const valence = own.reduce((sum, [, , type]) => sum + type, 0) + hhh - 1;
    if (valence > 14)
      continue;
    lines[atom.line] = line.padEnd(51).substring(0, 48) + String(valence).padStart(3) + line.substring(51);
    changed = true;
  }
  return changed ? lines.join('\n') : molBlock;
}

export function getQueryMolSafe(queryMolString: string, queryMolBlockFailover: string,
  rdKitModule: RDModule): RDMol | null {
  if (queryMolString && !hasNewLines(queryMolString) && queryMolString.length > MAX_SMILES_LENGTH)
    return null; // do not attempt to parse very long SMILES, will cause MOB.
  let queryMol = null;

  if (isMolBlock(queryMolString)) {
    // a query molblock is read as a query, unsanitized, V2000 or V3000: sanitizing changes what it asks for (an aromatic bond
    // outside an aromatic ring read single, crux-sketch spike query-roundtrip; an atom whose valence the molblock fixes read as
    // a radical, so that it matched radicals alone). One with an H count query on an atom (MDL's hhh, V3000's HCOUNT, which
    // Crux and Ketcher write) is read sanitized but not kekulized, unless RDKit cannot sanitize it: read unsanitized, its
    // Kekulé ring stayed Kekulé, so that it matched no aromatic ring (a V2000 one with its hhh atoms given their H, see
    // withHCountHs). So is one with an uncharged pentavalent N (a nitro drawn N(=O)=O, an azide N=N#N, an N-oxide N=O):
    // read unsanitized it matched nothing, the molecules holding those charge-separated, as sanitizing makes the query.
    // One with explicit H is read with its H merged into their atoms, which only the molecule reader does.
    if (queryMolString.includes(' H '))
      queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, rdKitModule).mol;
    else {
      if (hasHCountQuery(queryMolString) || hasPentavalentN(queryMolString)) {
        // with its hhh atoms given their H, else as it is (an H the atom has no room for: RDKit cannot sanitize it)
        for (const molBlock of new Set([withHCountHs(queryMolString), queryMolString])) {
          try {
            queryMol = rdKitModule.get_mol(molBlock, JSON.stringify({kekulize: false}));
          } catch (e) {
            queryMol = null;
          }
          if (queryMol !== null)
            return queryMol;
        }
      }
      try {
        queryMol = rdKitModule.get_qmol(queryMolString);
      } catch (e) {
        queryMol = null;
      }
      // Kekulé→aromatic conversion is a precaution for benzene-like queries; if
      // it throws on query-atom features (e.g. L + M ALS), keep the queryMol —
      // it's still usable for substructure matching, just won't aromatize.
      if (queryMol !== null) {
        try {
          queryMol.convert_to_aromatic_form();
        } catch (e) {/* non-fatal */}
      }
    }
  } else { // not a molblock
    const treatAsSmarts = _isSmarts(queryMolString);
    try {
      queryMol = rdKitModule.get_qmol(queryMolString);
    } catch (e) { }
    if (queryMol !== null) {
      const mol = rdKitModule.get_mol(queryMolString, '{"mergeQueryHs":true}');
      if (mol !== null) { // check the qmol is proper
        // for a plain SMILES use the molecule query: its pattern fingerprint and matching are consistent
        // with molblock-derived queries (a SMARTS qmol's pattern fp over-filters and drops valid hits).
        // keep the qmol for a real SMARTS unless it fails to match its own molecule, and whenever the SMILES
        // reading has radicals ([OH], [CH3], [F], [N+], c1cc[n+]cc1: meant as SMARTS)
        if (!hasRadicals(mol) && (!treatAsSmarts || mol.get_substruct_match(queryMol) === '{}')) {
          queryMol.delete(); //remove mol object previously stored in queryMol
          queryMol = mol;
        } else
          mol.delete();
      } // else, this looks to be a real SMARTS
    } else { // failover to queryMolBlockFailover
      // possibly get rid of fall-over in future
      queryMol = getMolSafe(queryMolBlockFailover, {mergeQueryHs: true}, rdKitModule).mol;
    }
  }
  return queryMol;
}

/**
 * The query a substructure search builds from `molString` ({@link getQueryMolSafe}), as RDKit's SMARTS: what names a query
 * where two drawings must be told apart exactly when the search tells them apart (the substructure filter's key). Null
 * when RDKit cannot read it.
 */
export function querySmarts(molString: string, rdKitModule: RDModule): string | null {
  const queryMol = getQueryMolSafe(molString, '', rdKitModule);
  try {
    return queryMol?.get_smarts() ?? null;
  } catch {
    return null;
  } finally {
    queryMol?.delete();
  }
}
