import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export enum MonomerWidthMode {
  long = 'long',
  short = 'short',
}

export const enum Tags {
  RendererSettingsChanged = '.mm.cellRenderer.settingsChanged',
}

export const rendererSettingsChangedState = {
  true: '1',
  false: '0',
};

export const enum Temps {
  monomerWidth = '.mm.cellRenderer.monomerWidth',
  colorCode = '.mm.cellRenderer.colorCode',
  compareWithCurrent = '.mm.cellRenderer.compareWithCurrent',
  highlightDifference = '.mm.cellRenderer.highlightDifference',
  gapLength = '.mm.cellRenderer.gapLength',
  monomerPlacer = '.mm.cellRenderer.monomerPlacer',
}

export const enum tempTAGS {
  referenceSequence = 'reference-sequence',
  currentWord = 'current-word',
  monomerWidthMode = 'monomer-width-mode',
}

// export const MacromoleculeCellRendererDefaults = new class {
//   monomerWidth: MonomerWidthMode = MonomerWidthMode.short;
//   maxMonomerLength: number = 3;
//   colorCode: boolean = true;
//   compareWithCurrent: boolean = true;
// }();
