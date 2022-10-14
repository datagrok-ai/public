import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {VdRegion} from '../vd-regions';
import {IViewer} from './viewer';

/** Interface for VdRegionsViewer from @datagrok/bio to unbind dependency to Bio package */
export interface IVdRegionsViewer extends IViewer {
  init(): Promise<void>;

  setDf(mlbDf: DG.DataFrame, regions: VdRegion[]): Promise<void>;
}