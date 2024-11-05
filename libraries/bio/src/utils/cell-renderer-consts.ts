import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const rendererSettingsChangedState = {
  true: '1',
  false: '0',
};

/** Macromolecule Cell Renderer Temp slot names */
export const enum MmcrTemps {
  maxMonomerLength = '.mm.cellRenderer.maxMonomerLength',
  colorCode = '.mm.cellRenderer.colorCode',
  compareWithCurrent = '.mm.cellRenderer.compareWithCurrent',
  highlightDifference = '.mm.cellRenderer.highlightDifference',
  gapLength = '.mm.cellRenderer.gapLength',
  monomerPlacer = '.mm.cellRenderer.monomerPlacer',

  rendererSettingsChanged = '.mm.cellRenderer.settingsChanged',
  overriddenLibrary = '.mm.cellRenderer.overriddenLibrary',
}

export const enum tempTAGS {
  referenceSequence = 'reference-sequence',
  currentWord = 'current-word',
}
