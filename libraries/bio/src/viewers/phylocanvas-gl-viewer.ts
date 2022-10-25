import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {IViewer} from './viewer';

export type NodeStyleType = { [propName: string]: any };
export type StylesType = { [nodeName: string]: NodeStyleType };

export interface IPhylocanvasGlViewer extends IViewer {
  get nwkDf(): DG.DataFrame;

  set nwkDf(value: DG.DataFrame);

  setProps(updater: { [propName: string]: any }): void;
}
