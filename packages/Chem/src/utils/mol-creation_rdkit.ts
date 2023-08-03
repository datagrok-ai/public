import {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {isMolBlock} from './chem-common';
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

export function isSmarts(molString: string): boolean {
  if (isMolBlock(molString))
    return MolfileHandler.getInstance(molString).isQuery();
  else
    return !!molString.match(/\[.?#\d|\$|&|;|,|!.?]/g);
}

export function getMolSafe(molString: string, details: object = {}, rdKitModule: RDModule,
  warnOff: boolean = true, checkIfSmarts: boolean = true): IMolContext {
  let isQMol = false;
  const kekulizeProp = (details as any).kekulize;
  let kekulize: boolean = typeof kekulizeProp === 'boolean' ? kekulizeProp : true;
  let useMolBlockWedging: boolean = false;
  let mol: RDMol | null = null;

  try {
    const _isSmarts = checkIfSmarts && isSmarts(molString);
    mol = _isSmarts ? rdKitModule.get_qmol(molString) : rdKitModule.get_mol(molString, JSON.stringify(details));
    isQMol = _isSmarts;
  } catch (e) {}
  if (!mol && kekulize) {
    kekulize = false; //Pyrrole cycles
    try {
      mol = rdKitModule.get_mol(molString, JSON.stringify({...details, kekulize}));
    } catch (e) {}
  }
  if (!mol) {
    isQMol = true;
    try {
      mol = rdKitModule.get_qmol(molString);
    } catch (e) {}
  }
  if (mol)
    useMolBlockWedging = (mol.has_coords() === 2);
  else if (!warnOff)
    console.error('Chem | In getMolSafe: RDKit.get_mol crashes on a molString: `' + molString + '`');
  return {mol, kekulize, isQMol, useMolBlockWedging};
}
