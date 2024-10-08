import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';
import {FilterSources, PositionHeight, WebLogoPropsDefault} from './web-logo';

// Data structures for V-Domain regions of antibodies

export enum VdRegionType {
  Unknown = 'unknown',
  FR = 'framework',
  CDR = 'cdr',
}

/** Describes V-DOMAIN (IG and TR) region (of multiple alignment)
 * https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
 * Attributes {@link positionStartName} and {@link positionEndName} are strings because
 * they correspond to position names as column names in ANARCI output (with a character index possible).
 */
export interface VdRegion {
  readonly type: VdRegionType;
  readonly name: string;
  readonly chain: string;

  readonly order: number;

  /** Name of the column containing sequence */
  readonly sequenceColumnName: string;

  /** Region start position (inclusive) */
  readonly positionStartName: string;

  /** Region end position (inclusive) */
  readonly positionEndName: string;
}

export const VdRegionsPropsDefault = new class {
  regionTypes: VdRegionType[] = [VdRegionType.CDR];
  chains: string[] = ['Heavy', 'Light'];
  //sequenceColumnNamePostfix: string = 'chain sequence';

  skipEmptyPositions: boolean = false;
  fitWidth: boolean = false;
  positionWidth: number = 16;
  positionHeight: PositionHeight = PositionHeight.Entropy;

  filterSource: FilterSources = WebLogoPropsDefault.filterSource;
}();

export type VdRegionsProps = Required<typeof VdRegionsPropsDefault>;

/** Interface for VdRegionsViewer from @datagrok/bio to unbind dependency to Bio package */
export interface IVdRegionsViewer extends VdRegionsProps, IViewer {
  init(): Promise<void>;

  /** List of regions referencing the column with sequences by {@link VdRegion.sequenceColumnName }. */
  setData(regions: VdRegion[]): void;
}

declare module 'datagrok-api/dg' {
  export interface DataFramePlotHelper {
    // eslint-disable-next-line max-len
    fromType(viewerType: 'VdRegions', options: Partial<VdRegionsProps>): Promise<DG.Viewer<VdRegionsProps> & IVdRegionsViewer>;
  }
}
