import {HelmType} from '../helm/types';

/** Type for user settings of monomer library set to use. */
export type UserLibSettings = {
  exclude: string[],
  explicit: string[],
  /** Libraries might contain different monomers for same symbol. Object for monomer symbol to library name of choice*/
  duplicateMonomerPreferences: {[polymerType: string]: {[monomerSymbol: string]: string}}
}

export type MonomerMapValue = { biotype: HelmType, symbol: string, atoms: number[], bonds: number[] };

export class MonomerMap extends Map<number, MonomerMapValue> {
  constructor(entries?: [number, MonomerMapValue][] | null) {
    super(entries);
  }
}

/** @property monomers key - the monomer's position in the sequence (gaps counted; HELM RNA counts sugar, base and
 * linker; a FASTA/separator nucleotide is one position) */
export class MolfileWithMap {
  constructor(
    public readonly molfile: string,
    public readonly monomers: MonomerMap,
  ) {}

  static createEmpty() { return new MolfileWithMap('', new MonomerMap(null)); }
}
