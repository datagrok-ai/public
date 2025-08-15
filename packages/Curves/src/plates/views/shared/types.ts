import * as DG from 'datagrok-api/dg';
import {Plate} from '../../../plate/plate';

export type PlateFile = {
  plate: Plate;
  file: DG.FileInfo;
  reconciliationMap: Map<string, string>;
  commonProperties?: Map<string, any>; // <-- ADD THIS PROPERTY
};

export type TemplateState = {
  plates: PlateFile[];
  activePlateIdx: number;
};

// ... rest of the file is unchanged
export interface PlateStateChangeEvent {
 type: 'plate-added' | 'plate-removed' | 'plate-selected' | 'mapping-changed' | 'template-changed';
 templateId?: number;
 plateIndex?: number;
 plate?: PlateFile;
}

export interface ValidationResult {
 element: HTMLElement;
 conflictCount: number;
}

export interface MappingDialogOptions {
 allPlates: PlateFile[];
 sourceMappings: Map<string, string>;
 onSync: (mappings: Map<string, string>, selectedIndexes: number[]) => void;
 onUndo: (mappedField: string) => void;
}
