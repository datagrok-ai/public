/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from '../../../plate/plate';

export type PlateFile = {
  plate: Plate;
  file: DG.FileInfo;
  commonProperties?: Map<string, any>;
};


export type TemplateState = {
  plates: PlateFile[];
  activePlateIdx: number;
};

export interface PlateStateChangeEvent {
 type: 'plate-added' | 'plate-removed' | 'plate-selected' | 'mapping-changed' | 'template-changed' | 'analysis-mapping-changed'| 'identifier-changed' | 'plate-data-changed'; // <-- ADD THIS
 templateId?: number;
 plateIndex?: number;
 plate?: PlateFile;
  analysisType?: 'drc' | 'doseRatio' | 'qpcr';
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
