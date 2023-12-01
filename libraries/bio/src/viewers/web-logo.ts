import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Observable} from 'rxjs';

import {IViewer} from './viewer';

export enum TAGS {
  /** Controls displaying WebLogo in a Macromolecule column's header tooltip */
  tooltipWebLogo = '.tooltipWebLogo',
}

export enum PositionHeight {
  Entropy = 'Entropy',
  full = '100%',
}

/** top, middle, bottom */
export enum VerticalAlignments {
  TOP = 'top',
  MIDDLE = 'middle',
  BOTTOM = 'bottom',
}

/** left, center, right */
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
  /** Sequences in selection are considered to plot WebLogo for faster exploration.
   * In case selection is empty displays all.
   */
  Selected = 'Selected',
}

export const WebLogoPropsDefault = new class {
  // -- Data --
  sequenceColumnName: string | null = null;
  /** Aggregation function for values of {@link valueColumnName} */
  valueAggrType: DG.AggregationType = DG.AGG.TOTAL_COUNT;
  /** Column name for values */
  valueColumnName: string = 'Activity';
  startPositionName: string | null = null;
  endPositionName: string | null = null;
  skipEmptySequences: boolean = true;
  skipEmptyPositions: boolean = false;
  shrinkEmptyTail: boolean = true;

  // -- Style --
  backgroundColor: number = 0xFFFFFFFF;
  positionHeight: string = PositionHeight.Entropy; // that is the way in the bioinformatics domain
  positionWidth: number = 16;

  // -- Layout --
  verticalAlignment: VerticalAlignments = VerticalAlignments.MIDDLE;
  horizontalAlignment: HorizontalAlignments = HorizontalAlignments.CENTER;
  fixWidth: boolean = false;
  fitArea: boolean = true;
  minHeight: number = 50;
  maxHeight: number = 100;
  showPositionLabels: boolean = true;
  positionMarginState: PositionMarginStates = PositionMarginStates.AUTO;
  positionMargin: number = 0;

  // -- Behavior --
  filterSource: FilterSources = FilterSources.Filtered;
}();

export type WebLogoProps = typeof WebLogoPropsDefault;

export interface IWebLogoViewer extends WebLogoProps, IViewer {
  get onSizeChanged(): Observable<void>;

  get positionMarginValue(): number;

  setOptions(options: Partial<WebLogoProps>): void;
}

export const positionRe: RegExp = /(\d+)([A-Z]?)/;

declare module 'datagrok-api/dg' {
  export interface DataFramePlotHelper {
    fromType(viewerType: 'WebLogo', options: Partial<WebLogoProps>): Promise<DG.Viewer<WebLogoProps> & IWebLogoViewer>;
  }
}
