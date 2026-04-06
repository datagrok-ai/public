import * as DG from 'datagrok-api/dg';

import {TwbDatasource, TwbColumn} from './tableau-types';


function toDgColumnType(datatype: string): DG.COLUMN_TYPE {
  switch (datatype) {
    case 'integer': return DG.COLUMN_TYPE.INT;
    case 'real': return DG.COLUMN_TYPE.FLOAT;
    case 'boolean': return DG.COLUMN_TYPE.BOOL;
    case 'date': return DG.COLUMN_TYPE.DATE_TIME;
    default: return DG.COLUMN_TYPE.STRING;
  }
}


function applyColumnTags(dgCol: DG.Column, col: TwbColumn): void {
  if (col.caption && col.caption !== col.name)
    dgCol.setTag('tableau.caption', col.caption);
  if (col.datatype)
    dgCol.setTag('tableau.datatype', col.datatype);
  if (col.role)
    dgCol.setTag('tableau.role', col.role);
  if (col.type)
    dgCol.setTag('tableau.type', col.type);
  if (col.aggregation)
    dgCol.setTag('tableau.aggregation', col.aggregation);
  if (col.containsNull)
    dgCol.setTag('tableau.containsNull', 'true');
  if (col.caption && col.caption !== col.name)
    dgCol.setTag('description', col.caption);
}


export function twbDatasourceToDataFrame(ds: TwbDatasource): DG.DataFrame {
  const df = DG.DataFrame.create(0);
  df.name = ds.caption || ds.name;

  df.setTag('source.format', 'tableau-twb');
  if (ds.sourceFile)
    df.setTag('tableau.sourceFile', ds.sourceFile);
  if (ds.sourceDirectory)
    df.setTag('tableau.sourceDirectory', ds.sourceDirectory);

  for (const col of ds.columns) {
    const dgCol = df.columns.addNew(col.name, toDgColumnType(col.datatype));
    applyColumnTags(dgCol, col);
  }

  return df;
}
