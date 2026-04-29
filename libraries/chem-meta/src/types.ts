import {ChemTemps} from './consts';


/** Persistent column tag used to link a SMILES (Molecule) column to a
 *  specific Molecule3D column for the Chem atom-picker bridge. The tag
 *  is set on the SMILES column; its value is the linked Molecule3D column
 *  name. Canonical constant owned by Chem's `constants.ts`; re-exported
 *  from the bio library so BiostructureViewer, Docking and other packages
 *  that already depend on `@datagrok-libraries/bio` can reference it
 *  without introducing a hard dependency on `@datagrok/chem`. Keep the
 *  string value in sync with Chem's `CHEM_ATOM_PICKER_LINKED_COL`. */
export const CHEM_ATOM_PICKER_LINKED_3D_COL_TAG = '.%chem-atom-picker-linked-3d-col';
// reverse of previous, this tag is set on ligand column
export const CHEM_ATOM_PICKER_LINKED_SMILES_COL = '.%chem-atom-picker-linked-smiles-col';

export const CHEM_ATOM_SELECTION_EVENT = 'd4-chem-atom-selection-changed';
//  viewer on 3D atom hover; Chem's renderer listens (reverse bridge). */
export const CHEM_MOL3D_SELECTION_EVENT = 'd4-chem-mol3d-selection-changed';

export interface Mol3DHoverEventArgs {
  mol3DColumnName: string;
  rowIdx: number;
  atom3DSerial: number | null;
  mode: 'preview' | 'paint' | 'erase';
}

// -- Types -------------------------------------------------------------------

/** Mirrors AtomIndexMapping from atom-index-mapper.ts (no cross-package import). */
export interface AtomMapping3D {
  mapping: number[];
  method: string;
  mappedCount: number;
  pdbSerials?: number[];
}

/** Shape of the cross-package `chem-interactive-selection-changed` event. */
export interface ChemSelectionEventArgs {
  column?: unknown;
  rowIdx: number;
  atoms: number[];
  mapping3D?: AtomMapping3D | null;
  persistent?: boolean;
  clearAll?: boolean;
  mol3DColumnName?: string;
}

export interface ISubstruct {
  atoms?: number[],
  bonds?: number[],
  highlightAtomColors?: { [key: number]: number[] | null },
  highlightBondColors?: { [key: number]: number[] | null },
  alignByScaffold?: string,
}

export function mergeSubstructs(substructs: ISubstruct[]): ISubstruct {
  const res: ISubstruct = {atoms: [], bonds: [], highlightAtomColors: {}, highlightBondColors: {}};
  for (const s of substructs) {
    res.atoms = [...res.atoms ?? [], ...s.atoms ?? []];
    res.bonds = [...res.bonds ?? [], ...s.bonds ?? []];
    res.highlightAtomColors = {...res.highlightAtomColors, ...s.highlightAtomColors};
    res.highlightBondColors = {...res.highlightBondColors, ...s.highlightBondColors};
  }
  return res;
}

export interface ISubstructProvider {
  /** To highlight */
  getSubstruct(tableRowIndex: number | null): ISubstruct | undefined;
}

/** Atom-picker variant of `ISubstructProvider` written by Chem's
 *  `AtomPickerController` into `col.temp[ChemTemps.SUBSTRUCT_PROVIDERS]` and
 *  read by BiostructureViewer's 3D Molstar viewer to mirror the highlighted
 *  atoms in 3D. The `__atomPicker` marker distinguishes these from other
 *  providers (MMP, scaffold, monomer-hover, etc.) that share the same
 *  `col.temp` array. Lives in chem-meta as the cross-package contract — Chem
 *  is the writer, BSV is the reader. */
export interface AtomPickerProvider extends ISubstructProvider {
  __atomPicker?: boolean;
  __rowIdx?: number;
  __atoms?: Set<number>;
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
