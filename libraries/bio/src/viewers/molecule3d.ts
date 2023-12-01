import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export type Molecule3DData = {
  /** Base64 encoded structure data (for text and binary) */ data: string;
  /** File extension, data format */ ext: string;
};

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
