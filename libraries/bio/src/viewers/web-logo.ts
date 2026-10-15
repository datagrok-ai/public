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

/** Functions with this role supply WebLogo profiles: `(table: dataframe, column: string, profile: string) => object`
 * returning a {@link WebLogoProfile}; `column` is the viewer's sequence column (e.g. the HMM package's
 * profile HMM logos). */
export const WEB_LOGO_PROFILE_ROLE = 'webLogoProfile';

/** A sequence profile drawn instead of the column's composition: per position,
 * the height of each monomer (e.g. a profile HMM's emission information in bits). */
export interface WebLogoProfile {
  /** Shown in tooltips, e.g. 'V-set (PF07686.20)' */
  title: string;
  /** Monomer colors: {@link ALPHABET} PT, DNA or RNA */
  alphabet: string;
  /** Units of the heights, e.g. 'bits' */
  units: string;
  positions: WebLogoProfilePosition[];
  /** For each position, the 0-based position of the sequence column aligned to it (-1: none), when the column
   * is aligned to the profile; clicking a letter then selects the rows that have it there. */
  columnPositions?: number[];
}

export interface WebLogoProfilePosition {
  name: string;
  /** Height of each monomer (non-negative) */
  heights: { [monomer: string]: number };
  /** Further values shown in the position tooltip */
  info?: { [name: string]: string };
}

export const WebLogoPropsDefault = new class {
  // -- Data --
  sequenceColumnName: string | null = null;
  /** Profile drawn instead of the column's composition, from a {@link WEB_LOGO_PROFILE_ROLE} function */
  profile: string | null = null;
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
  positionWidth: number = 20;

  // -- Layout --
  verticalAlignment: VerticalAlignments = VerticalAlignments.MIDDLE;
  horizontalAlignment: HorizontalAlignments = HorizontalAlignments.CENTER;
  fixWidth: boolean = false;
  fitArea: boolean = true;
  minHeight: number = 25;
  maxHeight: number | null = 300;
  maxMonomerLetters: number = 5;
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
