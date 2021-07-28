import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import { TREATMENT_ARM } from "../constants";

export function getUniqueValues(df: DG.DataFrame, colName: string) {
    const uniqueIds = new Set();
    let column = df.columns.byName(colName);
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++){
        let value = column.get(i);
        if(value && !column.isNone(i))
        uniqueIds.add(column.get(i));
    }
    return uniqueIds;
  }

  export function changeEmptyStringsToUnknown(df: DG.DataFrame, colName: string) {
    let column = df.columns.byName(colName);
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++){
        if(column.isNone(i)){
            column.set(i, 'Unknown');
        }
    }
  }


export function addTreatmentArm(df: DG.DataFrame, dm: DG.DataFrame, columnsToExtract: string[]) {
    let withArm = grok.data.joinTables(df, dm, [ 'USUBJID' ], [ 'USUBJID' ], columnsToExtract, [ TREATMENT_ARM ], DG.JOIN_TYPE.LEFT, false);
    changeEmptyStringsToUnknown(withArm, TREATMENT_ARM);
    return withArm;
}