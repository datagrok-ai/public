import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';


export interface IPhylocanvasGlViewer extends IViewer {
  get nwkDf(): DG.DataFrame;

  set nwkDf(value: DG.DataFrame);
}