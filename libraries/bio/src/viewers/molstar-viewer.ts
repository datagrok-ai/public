import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {IViewer} from './viewer';
import {PdbResDataFrameType} from '../pdb/pdb-helper';

export interface IMolstarViewer extends IViewer {
  setData(df: PdbResDataFrameType | DG.DataFrame): void;
}

export interface ISaguaroViewer extends IViewer {
  setData(df: PdbResDataFrameType | DG.DataFrame): void;
}
