import * as DG from 'datagrok-api/dg';
import {MpxProject} from './mpx-types';
import {mwxWorksheetToDataFrame} from '../mwx/mwx-to-dataframe';


/** Converts all worksheets in an MPX project to Datagrok DataFrames. */
export function mpxProjectToDataFrames(project: MpxProject): DG.DataFrame[] {
  return project.worksheets.map((ws) => mwxWorksheetToDataFrame(ws));
}
