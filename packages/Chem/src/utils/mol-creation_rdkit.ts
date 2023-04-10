import {RDModule, RDMol} from "@datagrok-libraries/chem-meta/src/rdkit-api";
import { isMolBlock } from "./chem-common";
import { MolfileHandler } from "@datagrok-libraries/chem-meta/src/parsing-utils/molfile-handler";
import { elementsTable } from "../constants";

export interface IMolContext {
  mol: RDMol | null; // null when molString is invalid
  kekulize: boolean;
  isQMol: boolean;
  useMolBlockWedging: boolean;
}

const molfileHandler = MolfileHandler.createInstance(`
Datagrok empty molecule

0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`);

export function isSmarts(molString: string): boolean {
  if (isMolBlock(molString)) {
    molfileHandler.init(molString);
    for (const atom of molfileHandler.atomTypes) {
      if (elementsTable.indexOf(atom) === -1) {
        return true;
      }
    }
    return false;
  } else
    return !!molString.match(/\[.?#\d|\$|&|;|,|!|:|\*.?\]/g);
}

export function getMolSafe(molString: string, details: object = {}, rdKitModule: RDModule, warnOff: boolean = true): IMolContext {
  let isQMol = false;
  let kekulize: boolean = true;
  let useMolBlockWedging: boolean = false;
  let mol: RDMol | null = null;

  try {
    const _isSmarts = isSmarts(molString);
    mol = _isSmarts ? rdKitModule.get_qmol(molString) : rdKitModule.get_mol(molString, JSON.stringify(details));
    isQMol = _isSmarts;
  }
  catch (e) {
    if (mol !== null) {
      mol.delete();
      mol = null;
    }

    kekulize = false;
    try { mol = rdKitModule.get_mol(molString, JSON.stringify({ ...details, kekulize })); }
    catch (e2) {
      if (mol !== null) {
        mol.delete();
        mol = null;
      }

      try { mol = rdKitModule.get_qmol(molString); }
      catch (e3) {
        if (mol !== null) {
          mol.delete();
          mol = null;
        }
        if (!warnOff)
          console.error('Chem | In getMolSafe: RDKit.get_mol crashes on a molString: `' + molString + '`');
        return { mol, kekulize, isQMol, useMolBlockWedging };
      }
      return { mol, kekulize, isQMol, useMolBlockWedging };
    }
    if (mol.is_valid()) {
      useMolBlockWedging = mol.has_coords();
    }
    return { mol, kekulize, isQMol, useMolBlockWedging };
  }
  if (mol.is_valid()) {
    useMolBlockWedging = mol.has_coords();
  } else {
    mol?.delete()
    mol = null;
  }
  return { mol, kekulize, isQMol, useMolBlockWedging };
}
