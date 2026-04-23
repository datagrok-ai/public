import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export type Molecule3DData = {
  /** Base64 encoded structure data (for text and binary) */ data: string;
  /** File extension, data format */ ext: string;
};

/** Persistent column tag used to link a SMILES (Molecule) column to a
 *  specific Molecule3D column for the Chem atom-picker bridge. The tag
 *  is set on the SMILES column; its value is the linked Molecule3D column
 *  name. Canonical constant owned by Chem's `constants.ts`; re-exported
 *  from the bio library so BiostructureViewer, Docking and other packages
 *  that already depend on `@datagrok-libraries/bio` can reference it
 *  without introducing a hard dependency on `@datagrok/chem`. Keep the
 *  string value in sync with Chem's `CHEM_ATOM_PICKER_LINKED_COL`. */
export const CHEM_ATOM_PICKER_LINKED_COL = '%chem-atom-picker-linked-col';

export enum DockingTags {
  dockingRole = 'docking.role',
  dockingTarget = 'docking.target',
}

export enum DockingRole {
  target = 'target',
  ligand = 'ligand',
};

export type DockingTargetData = Molecule3DData;
export type DockingLigandData = Molecule3DData;

export interface IMolecule3DBrowser {

  showStructure(data: Molecule3DData): void;
}
