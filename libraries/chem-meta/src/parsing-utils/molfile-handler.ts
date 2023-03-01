/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AbstractChemicalTableParser, AtomAndBondCounts} from './chemical-table-parser';

enum MOLFILE_VERSION {
  V2000 = 'V2000',
  V3000 = 'V3000',
}

// todo: port exportable constants to a separate module
const enum V3K {
  HEADER = '999 V3000\n',
  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  BEGIN_ATOM_BLOCK = 'M  V30 BEGIN ATOM',
  END = 'M  END',
}

const enum V2K {
  HEADER = '999 V2000\n',
  END = 'M  END',
}

/** Class for handling MDL Molfiles V2000 and V3000 */
class MolfileHandler extends AbstractChemicalTableParser {
  constructor(molfile: string) {
    super(molfile);
    this.molfileVersion = MolfileHandler.determineMolfileVersion(molfile);
  }

  private molfileVersion: MOLFILE_VERSION;

  /** Determine whether the file is V2000/V3000, or throw */
  private static determineMolfileVersion(molfile: string): MOLFILE_VERSION {
    if (MolfileHandler.validateV3K(molfile))
      return MOLFILE_VERSION.V3000;
    else if (MolfileHandler.validateV2K(molfile))
      return MOLFILE_VERSION.V3000;
    else
      throw new Error('Malformed molfile');
  }

  protected parseAtomAndBondCounts(): AtomAndBondCounts {
    this.currentIdx = this.getAtomBlockIdx();
    if (this.molfileVersion === MOLFILE_VERSION.V2000)
      this.parseAtomAndBondCountsV2K();
    else
      this.parseAtomAndBondCountsV3K();
  };

  protected parseAtomCoordinates(): Float32Array[] {

  };

  protected parseAtomTypes(): string[] {

  };

  protected getAtomBlockIdx(): number {
    if (this.molfileVersion === MOLFILE_VERSION.V2000)
      return this.getAtomBlockIdxV2K();
    else
      return this.getAtomBlockIdxV3K();
  };

  private getAtomBlockIdxV2K(): number {
    let idx = 0;
    // 4 is for 3 header lines in V2K + one counts line
    for (let i = 0; i < 4; ++i)
      idx = this.getIdxOfNextLine(idx);
    return idx;
  }

  private getAtomBlockIdxV3K(): number {
    let idx = this.file.indexOf(V3K.BEGIN_ATOM_BLOCK);
    idx = this.getIdxOfNextLine(idx);
    return idx;
  }

  protected getBondBlockIdx(): number {
    if (this.molfileVersion === MOLFILE_VERSION.V2000)
      return this.getBondBlockIdxV2K();
    else
      return this.getBondBlockIdxV3K();
  };

  private getBondBlockIdxV2K(): number {
    let idx = this.getAtomBlockIdx();
    for (let i = 0; i < )
  }

  // todo: devise a more complex validation
  private static validateV3K(molfile: string): boolean {
    if (molfile.indexOf(V3K.HEADER) !== -1 &&
    molfile.indexOf(V3K.END) !== -1)
      return true;
    return false;
  }

  // todo: devise a more complex validation
  private static validateV2K(molfile: string): boolean {
    if (molfile.indexOf(V2K.HEADER) !== -1 &&
    molfile.indexOf(V2K.END) !== -1)
      return true;
    return false;
  }

  // private setAtomAndBondCounts(): void {
  //   if (!this.atomCount) {
  //     const parse = (this.molfileVersion === MolfileVersion.V2000) ? this.parseAtomAndBondCountsV3K : this.parseAtomAndBondCountsV2K;
  //     const counts = parse();
  //     this.atomCount = counts.atomCount;
  //     this.bondCount = counts.bondCount;
  //   }
  // }

  private parseAtomAndBondCountsV2K(): AtomAndBondCounts {

    // parse atom count
    let begin = this.molfile.indexOf(V3K_BEGIN_COUNTS_LINE) + V3K_COUNTS_SHIFT;
    let end = molfileV3K.indexOf(' ', begin + 1);
    const numOfAtoms = parseInt(molfileV3K.substring(begin, end));

    // parse bond count
    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin + 1);
    const numOfBonds = parseInt(molfileV3K.substring(begin, end));

    return {atomCount: numOfAtoms, bondCount: numOfBonds};
  };

  private parseAtomAndBondCountsV3K(): AtomAndBondCounts {
    molfileV3K = molfileV3K.replaceAll('\r', ''); // to handle old and new sdf standards

    // parse atom count
    let begin = molfileV3K.indexOf(V3K_BEGIN_COUNTS_LINE) + V3K_COUNTS_SHIFT;
    let end = molfileV3K.indexOf(' ', begin + 1);
    const numOfAtoms = parseInt(molfileV3K.substring(begin, end));

    // parse bond count
    begin = end + 1;
    end = molfileV3K.indexOf(' ', begin + 1);
    const numOfBonds = parseInt(molfileV3K.substring(begin, end));

    return {atomCount: numOfAtoms, bondCount: numOfBonds};
  }
}
