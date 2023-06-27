import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export enum MonomerWidthMode {
  long = 'long',
  short = 'short',
}

export const enum Tags {
  calculated = '.mm.cellRenderer.calculated',
}

export const enum Temps {
  monomerWidth = '.mm.cellRenderer.monomerWidth',
  maxMonomerLength = '.mm.cellRenderer.maxMonomerLength',
  colorCode = '.mm.cellRenderer.colorCode',
  compareWithCurrent = '.mm.cellRenderer.compareWithCurrent',
  highlightDifference = '.mm.cellRenderer.highlightDifference',
}

// export const MacromoleculeCellRendererDefaults = new class {
//   monomerWidth: MonomerWidthMode = MonomerWidthMode.short;
//   maxMonomerLength: number = 3;
//   colorCode: boolean = true;
//   compareWithCurrent: boolean = true;
// }();
