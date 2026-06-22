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

export function _isSmarts(molString: string): boolean {
  if (isMolBlock(molString))
    return MolfileHandler.getInstance(molString).isQuery();
  else
    return !!molString.match(/\[.?#\d|\$|&|;|,|!.?]/g);
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


export function getQueryMolSafe(queryMolString: string, queryMolBlockFailover: string,
  rdKitModule: RDModule): RDMol | null {
  if (queryMolString && !hasNewLines(queryMolString) && queryMolString.length > MAX_SMILES_LENGTH)
    return null; // do not attempt to parse very long SMILES, will cause MOB.
  let queryMol = null;

  if (isMolBlock(queryMolString)) {
    if (queryMolString.includes(' H ') || queryMolString.includes('V3000'))
      queryMol = getMolSafe(queryMolString, {mergeQueryHs: true}, rdKitModule).mol;
    else {
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
        // keep the qmol only for a real SMARTS, or when the qmol fails to match its own molecule
        if (!treatAsSmarts || mol.get_substruct_match(queryMol) === '{}') {
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
