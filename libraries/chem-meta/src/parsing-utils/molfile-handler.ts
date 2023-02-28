/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

enum MolfileVersion {
  V2000 = 'V2000',
  V3000 = 'V3000',
}

const enum V3K {
  HEADER = '999 V3000\n',
  BEGIN_COUNTS_LINE = 'M  V30 COUNTS ',
  END = 'M  END',
}

const enum V2K {
  HEADER = '999 V2000\n',
  END = 'M  END',
}

type Coordinates = {
  x: Float32Array,
  y: Float32Array,
  z: Float32Array,
}

/** Class for handling MDL Molfiles V2000 and V3000 */
class MolfileHandler {
  constructor(molfile: string) {
    this.molfileVersion = MolfileHandler.determineMolfileVersion(molfile);
    this.molfile = molfile.replaceAll('\r', ''); // to handle old and new sdf standards
    this.atomCoordinates = null;
    this.atomCount = null;
    this.bondCount = null;
    this.atomTypes = null;
    this.currentIdx = 0;
  }

  private molfile: string;
  private molfileVersion: MolfileVersion;
  private atomCount: null | number;
  private bondCount: null | number;
  /** X, Y and Z coordinates of atoms  */
  private atomCoordinates: null | Coordinates;
  private atomTypes: null | string[];
  private currentIdx: number;

  /** Determine whether the file is V2000/V3000, or throw */
  private static determineMolfileVersion(molfile: string): MolfileVersion {
    if (MolfileHandler.validateV3K(molfile))
      return MolfileVersion.V3000;
    else if (MolfileHandler.validateV2K(molfile))
      return MolfileVersion.V3000;
    else
      throw new Error('Malformed molfile');
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

  private parseCoordinates(): Float32Array[] {
    const x = new Float32Array();
    const y = new Float32Array();
    const z = new Float32Array();
    return [x, y, z];
  }

  private setAtomAndBondCounts(): void {
    if (!this.atomCount) {
      const parse = (this.molfileVersion === MolfileVersion.V2000) ? this.parseAtomAndBondCountsV3K : this.parseAtomAndBondCountsV2K;
      const counts = parse();
      this.atomCount = counts.atomCount;
      this.bondCount = counts.bondCount;
    }
  }

  private jumpToNextLine(): void {
    this._currentIdx = this.getIdxOfNextLine();
  }

  private getIdxOfNextLine(): number {
    if (this._str.at(this._currentIdx) !== '\n')
      return this._str.indexOf('\n', this._currentIdx) + 1;
    else
      return this._str.indexOf('\n', this._currentIdx + 1) + 1;
  }

  private parseAtomAndBondCountsV2K() { atomCount: number, bondCount: number } {

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

  private parseAtomAndBondCountsV3K() { atomCount: number, bondCount: number } {
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

  public get x: Float32Array {
    if (!this.) {

    } else {
      return this.x;
    }
  }

}
