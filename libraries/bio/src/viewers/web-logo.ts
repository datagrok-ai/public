import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';

export enum PositionHeight {
  Entropy = 'Entropy',
  full = '100%',
}

export enum VerticalAlignments {
  TOP = 'top',
  MIDDLE = 'middle',
  BOTTOM = 'bottom',
}

export enum HorizontalAlignments {
  LEFT = 'left',
  CENTER = 'center',
  RIGHT = 'right',
}

export enum PositionMarginStates {
  AUTO = 'auto',
  ON = 'on',
  OFF = 'off',
}

/** The source for filtering sequences considered to plot WebLogo  */
export enum FilterSources {
  /** Sequences of filtered rows are considered, default. */
  Filtered = 'Filtered',
  /** Sequences in selection are considered to plot WebLogo for faster exploration. */
  Selected = 'Selected',
}

export const WebLogoPropsDefault = new class {
  // -- Data --
  sequenceColumnName: string | null = null;
  startPositionName: string | null = null;
  endPositionName: string | null = null;
  skipEmptySequences: boolean = true;
  skipEmptyPositions: boolean = false;
  shrinkEmptyTail: boolean = true;

  // -- Style --
  backgroundColor: number = 0xFFFFFFFF;
  positionHeight: string = PositionHeight.full;
  positionWidth: number = 16;

  // -- Layout --
  verticalAlignment: VerticalAlignments = VerticalAlignments.MIDDLE;
  horizontalAlignment: HorizontalAlignments = HorizontalAlignments.CENTER;
  fixWidth: boolean = false;
  fitArea: boolean = true;
  minHeight: number = 50;
  maxHeight: number = 100;
  positionMarginState: PositionMarginStates = PositionMarginStates.AUTO;
  positionMargin: number = 0;

  // -- Behavior --
  filterSource: FilterSources = FilterSources.Filtered;
}();

export type WebLogoProps = Required<typeof WebLogoPropsDefault>;

export interface IWebLogoViewer extends WebLogoProps, IViewer {
  // Some methods to control the viewer
}

export const positionSeparator: string = ', ';
export const positionRe: RegExp = /(\d+)([A-Z]?)/;

export enum TAGS {
  positionNames = '.positionNames',
}

declare module 'datagrok-api/dg' {
  export interface DataFramePlotHelper {
    fromType(viewerType: 'WebLogo', options: Partial<WebLogoProps>): Promise<DG.Viewer & IWebLogoViewer>;
  }
}
