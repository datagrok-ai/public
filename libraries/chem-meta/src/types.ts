import {ChemTemps} from './consts';

export interface ISubstruct {
  atoms?: number[],
  bonds?: number[],
  highlightAtomColors?: { [key: number]: number[] | null },
  highlightBondColors?: { [key: number]: number[] | null },
  alignByScaffold?: string,
  /** Per-atom replacement labels (e.g. {0: 'X'}). Forwarded to RDKit MinimalLib
   *  via the `atomLabels` field in the draw details JSON. Replaces the rendered
   *  atom symbol for the given atom indices. */
  atomLabels?: { [key: number]: string },
  /** Per-atom annotation strings rendered next to atoms (e.g. pKa, partial
   *  charge, similarity contribution). Picked up by SVG-overlay renderers and
   *  by RDKit MinimalLib builds that support `atomNotes` in the details JSON. */
  atomNotes?: { [key: number]: string },
  /** Per-bond annotation strings rendered next to bonds. Same semantics as
   *  atomNotes. */
  bondNotes?: { [key: number]: string },
}

export function mergeSubstructs(substructs: ISubstruct[]): ISubstruct {
  const res: ISubstruct = {
    atoms: [], bonds: [],
    highlightAtomColors: {}, highlightBondColors: {},
    atomLabels: {}, atomNotes: {}, bondNotes: {},
  };
  for (const s of substructs) {
    res.atoms = [...res.atoms ?? [], ...s.atoms ?? []];
    res.bonds = [...res.bonds ?? [], ...s.bonds ?? []];
    res.highlightAtomColors = {...res.highlightAtomColors, ...s.highlightAtomColors};
    res.highlightBondColors = {...res.highlightBondColors, ...s.highlightBondColors};
    res.atomLabels = {...res.atomLabels, ...s.atomLabels};
    res.atomNotes = {...res.atomNotes, ...s.atomNotes};
    res.bondNotes = {...res.bondNotes, ...s.bondNotes};
  }
  return res;
}

export interface ISubstructProvider {
  /** To highlight */
  getSubstruct(tableRowIndex: number | null): ISubstruct | undefined;
}

export type MonomerHoverData = {
  dataFrameId: string,
  gridRowIdx: number,
  seqColName: string,
  seqPosition: number
  gridCell: any | null,
  /** Contains color of the monomer, empty lists on monomer that does not exist in molecule. */
  getSubstruct(): ISubstruct | undefined;
}

type MonomerHoverWindow = Window & {
  $monomerHover: MonomerHoverData | null;
}

declare const window: MonomerHoverWindow;

/** Return global monomer hover object. null - no monomer hover, negative seqPosition - hovered on not in a*/
export function getMonomerHover(): MonomerHoverData | null {
  return window.$monomerHover ?? null;
}

export function setMonomerHover(value: MonomerHoverData | null): void {
  window.$monomerHover = value;
}

export function addSubstructProvider(colTemp: any, substructProvider: ISubstructProvider): void {
  let list = colTemp[ChemTemps.SUBSTRUCT_PROVIDERS];
  if (!list)
    list = colTemp[ChemTemps.SUBSTRUCT_PROVIDERS] = [];
  list.push(substructProvider);
  colTemp[ChemTemps.SUBSTRUCT_PROVIDERS] = list;
}

export function getSubstructProviders(colTemp: any): ISubstructProvider[] {
  return colTemp?.[ChemTemps.SUBSTRUCT_PROVIDERS] ?? [];
}
