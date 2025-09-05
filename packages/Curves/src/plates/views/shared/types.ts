/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import {Plate} from '../../../plate/plate';

export type PlateFile = {
  plate: Plate;
  file: DG.FileInfo;
  reconciliationMap: Map<string, string>;
  analysisMappings: {
    drc: Map<string, string>;
    doseRatio: Map<string, string>;
  };
  commonProperties?: Map<string, any>;
};


export type TemplateState = {
  plates: PlateFile[];
  activePlateIdx: number;
};

// Add new event types for analysis mapping changes
export interface PlateStateChangeEvent {
  type: 'plate-added' | 'plate-removed' | 'plate-selected' | 'mapping-changed' | 'template-changed' | 'analysis-mapping-changed';
  templateId?: number;
  plateIndex?: number;
  plate?: PlateFile;
  analysisType?: 'drc' | 'doseRatio'; // New field for analysis-specific events
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
