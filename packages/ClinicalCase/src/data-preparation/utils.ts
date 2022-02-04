import * as DG from "datagrok-api/dg";
import * as grok from 'datagrok-api/grok';
import { VISIT_DAY, VISIT_NAME, SUBJECT_ID, TREATMENT_ARM } from "../columns-constants";

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

  export function filterNulls(df: DG.DataFrame, colName: string) {
    let column = df.columns.byName(colName);
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++){
        if(column.isNone(i)){
            df.rows.removeAt(i);
            i--;
            rowCount-=1;
        }
    }
  }

  export function filterFloats(df: DG.DataFrame, colName: string) {
    let column = df.columns.byName(colName);
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++){
        if(parseFloat(column.get(i)) === NaN){
            df.rows.removeAt(i);
            i--;
            rowCount-=1;
        }
    }
  }

  export function filterBooleanColumn(df: DG.DataFrame, colName: string, value: boolean) {
    let column = df.columns.byName(colName);
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++){
        if(column.get(i) === value){
            df.rows.removeAt(i);
            i--;
            rowCount-=1;
        }
    }
  }


  export function dateDifferenceInDays(start: string, end: string) {
    const startDate = new Date(start) as any;
    const endDate = new Date(end) as any;
    const diffTime = endDate - startDate;
    return Math.ceil(diffTime / (1000 * 60 * 60 * 24)); 
  }


export function addDataFromDmDomain(df: DG.DataFrame, dm: DG.DataFrame, columnsToExtract: string[], columnsToExtractFromDm: string[], subjIdColName = SUBJECT_ID) {
    let withArm = grok.data.joinTables(df, dm, [ subjIdColName ], [ SUBJECT_ID ], columnsToExtract, columnsToExtractFromDm, DG.JOIN_TYPE.LEFT, false);
   // columnsToExtractFromDm.forEach(it => changeEmptyStringsToUnknown(withArm, it));
    return withArm;
}


export function createFilteredTable(df: DG.DataFrame, groupCols: string[], condition: string) {
    return df
        .groupBy(groupCols)
        .where(condition)
        .aggregate();
}

export function dataframeContentToRow(df: DG.DataFrame) {
    let content = '';
    let rowCount = df.rowCount;
    for (let i = 0; i < rowCount; i++) {
        for (let column of df.columns) {
            content = `${content}${df.get(column.name, i)} `
        }
        content = `${content}; `;
    }
    return content;
}


export function getNullOrValue(df: DG.DataFrame, colname: string, index: number){
    return df.getCol(colname).isNone(index) ?  'null' : df.get(colname, index);
}


export function dictToString(dict: any){
    let str = '';
    Object.keys(dict).forEach(key => str+=`${key}: ${dict[key]}<br/>`);
    return str;
}

export function getVisitNamesAndDays(df: DG.DataFrame, allowNulls = false) {
    const data = df
        .groupBy([VISIT_DAY, VISIT_NAME])
        .aggregate();
    const visitsArray = [];
    let rowCount = data.rowCount;
    for (let i = 0; i < rowCount; i++) {
        const day = data.getCol(VISIT_DAY).isNone(i) ? null : data.get(VISIT_DAY, i);
        const name = data.getCol(VISIT_NAME).isNone(i) ? null : data.get(VISIT_NAME, i);
        if (allowNulls) {
            if (!visitsArray.filter(it => it.day === day && it.name === name).length) {
                visitsArray.push({day: day, name: name})
            };

        } else {
            if (day && name) {
                visitsArray.push({ day: day, name: name });
            }
        }

    }
    //@ts-ignore
    return visitsArray.sort((a, b) => { return (b.day !== null) - (a.day !== null) || a.day - b.day; });
}


  export function createPivotedDataframe(df: DG.DataFrame, groupByCols: string[], pivotCol: string, aggregatedColName: string, splitBy: string[]) {
    return df
        .groupBy(groupByCols.concat(splitBy))
        .pivot(pivotCol)
        .first(aggregatedColName)
        .aggregate();
}

export function createPivotedDataframeAvg(df: DG.DataFrame, groupByCols: string[], pivotCol: string, aggregatedColName: string, splitBy: string[]) {
    return df
        .groupBy(groupByCols.concat(splitBy))
        .pivot(pivotCol)
        .avg(aggregatedColName)
        .aggregate();
}

export function checkDateFormat(colToCheck: DG.Column, rowCount: number) {
    const rowsWithIncorrectDates = [];
    for (let i = 0; i < rowCount; i++) {
        if (!colToCheck.isNone(i) && isNaN(Date.parse(colToCheck.get(i))))
            rowsWithIncorrectDates.push(i);
    }
    return rowsWithIncorrectDates;
}

export function convertColToString(df: DG.DataFrame, col: string) {
    if (df.col(col).type !== 'string') {
        df.columns.addNewString(`${col}_s`).init((i) => df.get(col, i).toString())
        df.columns.remove(col);
        df.col(`${col}_s`).name = col;
    }
}