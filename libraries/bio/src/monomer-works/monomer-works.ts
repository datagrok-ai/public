import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import '../types/helm';
import * as org from 'org';

import {IMonomerLib} from '../types';

/** Hypothetical interface to convert mol block notation.
 * It should be placed in the chem-meta package, and have an implementation in the Chem package.
 * So dependency of MonomerWorks on molfile conversion operation becomes explicit.
 */
export interface IMolfileConverter {
  convertV2000toV3000(src: string): string;

  convertV3000toV2000(src: string): string;
}

export class MonomerWorks {
  private monomerLib: IMonomerLib;

  //private molfileConverter: IMolfileConverter;

  constructor(monomerLib: IMonomerLib/*, molfileConverter: IMolfileConverter*/) {
    this.monomerLib = monomerLib;
    //this.molfileConverter = molfileConverter;
  }

  //types according to Monomer possible
  public getCappedRotatedMonomer(monomerType: org.helm.PolymerType, monomerName: string): string | null {
    // TODO: Check type of monomerType arg
    const monomer = this.monomerLib.getMonomer(monomerType, monomerName);
    if (monomer)
      return monomer.molfile; //TODO cap

    return null;
  }
}

export function helmTypeToPolymerType(helmType: org.helm.HelmType): org.helm.PolymerType {
  let polymerType: org.helm.PolymerType | undefined = undefined;
  switch (helmType) {
  case 'HELM_BASE':
  case 'HELM_SUGAR': // r - ribose, d - deoxyribose
  case 'HELM_LINKER': // p - phosphate
  case 'HELM_NUCLETIDE':
    polymerType = 'RNA';
    break;
  case 'HELM_AA':
    polymerType = 'PEPTIDE';
    break;
  case 'HELM_CHEM':
    polymerType = 'CHEM';
    break;
  case 'HELM_BLOB':
    polymerType = 'BLOB';
    break;
  default:
    polymerType = 'PEPTIDE';
    console.warn(`Unexpected HelmType '${helmType}'`);
  }
  return polymerType;
}

export function getRS(smiles: string) {
  const newS = smiles.match(/(?<=\[)[^\][]*(?=])/gm);
  const res: { [name: string]: string } = {};
  let el = '';
  let digit;
  if (!!newS) {
    for (let i = 0; i < newS.length; i++) {
      if (newS[i] != null) {
        if (/\d/.test(newS[i])) {
          digit = newS[i][newS[i].length - 1];
          newS[i] = newS[i].replace(/[0-9]/g, '');
          for (let j = 0; j < newS[i].length; j++) {
            if (newS[i][j] != ':')
              el += newS[i][j];
          }
          res['R' + digit] = el;
          el = '';
        }
      }
    }
  }
  return res;
}
