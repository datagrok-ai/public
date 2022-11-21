import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {sequenceToMolFileST} from './to-atomic-level';

import {IMonomerLib, Monomer} from '../types';


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
  public getCappedRotatedMonomer(monomerType: string, monomerName: string): string | null {
    const monomer = this.monomerLib.getMonomer(monomerType, monomerName);
    if (monomer)
      return monomer.molfile; //TODO cap

    return null;
  }

  /** Consumes a list of monomer symbols and restores molfileV3K (SequenceTranslator/ST version) */
  public getAtomicLevel(monomers: string[], polymerType: string): string | null {
    return sequenceToMolFileST(
      monomers, this.monomerLib.getMonomerMolsByType(polymerType)!
    );
  }
}
