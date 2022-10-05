import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {VdRegion} from '../vd-regions';

/** Interface for VdRegionsViewer from @datagrok/bio to unbind dependency to Bio package */
export interface IVdRegionsViewer {
  get root(): HTMLElement;

  init(): Promise<void>;

  setDf(mlbDf: DG.DataFrame, regions: VdRegion[]): Promise<void>;
}