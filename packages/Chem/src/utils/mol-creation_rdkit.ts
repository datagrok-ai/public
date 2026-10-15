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
  warnOff: boolean = true): IMolContext {
  if (molString && !hasNewLines(molString) && molString.length > MAX_SMILES_LENGTH)
    return {mol: null, kekulize: true, isQMol: false, useMolBlockWedging: false}; // do not attempt to parse very long SMILES, will cause MOB.
  let isQMol = false;
  const kekulizeProp = (details as any).kekulize;
  let kekulize: boolean = typeof kekulizeProp === 'boolean' ? kekulizeProp : true;
  let useMolBlockWedging: boolean = false;
  let mol: RDMol | null = null;

  try {
    mol = rdKitModule.get_mol(molString, JSON.stringify(details));
  } catch (e) {}
  if (!mol && kekulize) {
    kekulize = false; //Pyrrole cycles
    try {
      mol = rdKitModule.get_mol(molString, JSON.stringify({...details, kekulize}));
    } catch (e) {}
  }
  if (!mol) {
    kekulize = true; //previous step failed, so retrun kekulize parameter to true
    isQMol = true;
    try {
      mol = rdKitModule.get_qmol(molString);
      isQMol = true;
    } catch (e) {}
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

export function getQueryMolSafe(queryMolString: string, queryMolBlockFailover: string,
  rdKitModule: RDModule): RDMol | null {
  if (queryMolString && !hasNewLines(queryMolString) && queryMolString.length > MAX_SMILES_LENGTH)
    return null; // do not attempt to parse very long SMILES, will cause MOB.
  let queryMol = null;

  if (isMolBlock(queryMolString)) {
    // a query molblock is read as a query, unsanitized, V2000 or V3000: sanitizing changes what it asks for (an aromatic bond
    // outside an aromatic ring read single, crux-sketch spike query-roundtrip; an atom whose valence the molblock fixes read as
    // a radical, so that it matched radicals alone; a pentavalent N drawn read charge-separated). One with an H count query on
    // an atom (MDL's hhh, V3000's HCOUNT, which Crux and Ketcher write) is read sanitized but not kekulized, unless RDKit cannot
    // sanitize it: read unsanitized, its Kekulé ring stayed Kekulé, so that it matched no aromatic ring. One with explicit H is
    // read with its H merged into their atoms, which only the molecule reader does.
    if (queryMolString.includes(' H '))
      queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, rdKitModule).mol;
    else {
      if (hasHCountQuery(queryMolString)) {
        try {
          queryMol = rdKitModule.get_mol(queryMolString, JSON.stringify({kekulize: false}));
        } catch (e) {
          queryMol = null;
        }
        if (queryMol !== null)
          return queryMol;
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
