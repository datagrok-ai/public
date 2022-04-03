
import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import { createTotalValuesForRowsAndCols, replaceNullsWithValues } from '../data-preparation/utils';
var { jStat } = require('jstat')

export function getDegreesOfFreedom(df: DG.DataFrame, categoriesCol: string, groupCol: string) {
    return (df.col(categoriesCol).categories.length - 1) * (df.col(groupCol).categories.length - 1);
}

export function createContingencyTable(df: DG.DataFrame, categoriesColumn: string, groupCol: string) {
    const countDf = df.groupBy([categoriesColumn].concat([groupCol]))
        .count()
        .aggregate();

    const contingencyTable = countDf
        .groupBy([groupCol])
        .pivot(categoriesColumn)
        .first('count')
        .aggregate();

    replaceNullsWithValues(contingencyTable, 0)
    return contingencyTable;
}


export function createContingencyTableWithExpectedVals(contTable: DG.DataFrame, groupCol: string, totalCol: string) {
    const rowCount = contTable.rowCount;
    const totalEvents = contTable.get(totalCol, rowCount - 1);

    for (let i = 0; i < rowCount - 1; i++) {
        let totalEventsInRow = contTable.get(totalCol, i);
        for (let column of contTable.columns) {
            if (column.name !== groupCol && column.name !== totalCol) {
                const totalEventsInCol = contTable.get(column.name, rowCount - 1);
                const expected = totalEventsInRow * totalEventsInCol / totalEvents;
                contTable.set(column.name, i, expected);
            }
        }
    }

    return contTable;
}


export function pearsonChiCriterion(contTable: DG.DataFrame, contTableWithExp: DG.DataFrame, groupCol: string, totalCol: string) {
    const rowCount = contTable.rowCount;
    let chiCriterion = 0;
    for (let i = 0; i < rowCount - 1; i++) {
        for (let column of contTable.columns) {
            if (column.name !== groupCol && column.name !== totalCol) {
                const actualVal = contTable.get(column.name, i);
                const expectedVal = contTableWithExp.get(column.name, i);
                chiCriterion += Math.pow((actualVal - expectedVal), 2) / expectedVal;
            }
        }
    }
    return chiCriterion;
}

export function getPearsonChiCriterionPValue(df: DG.DataFrame, categoryCol: string, groupCol: string){
    const totalColName = 'total';
    const contTable = createContingencyTable(df, categoryCol, groupCol);
    const contTableWithTotal = createTotalValuesForRowsAndCols(contTable, groupCol, totalColName);    
    const contTableWithExp = createContingencyTableWithExpectedVals(contTableWithTotal, groupCol, totalColName);   
    const chiCrit = pearsonChiCriterion(contTable, contTableWithExp, 'SEX', totalColName);    
    const degreesOfFreedom = getDegreesOfFreedom(df, categoryCol, groupCol);
    console.log(chiCrit);
    console.log(degreesOfFreedom);
    const pValue = jStat.chisquare.cdf( chiCrit, degreesOfFreedom );
    return pValue;
}
