import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { PERSON_ID } from '../constants';
import { cohorts } from '../cohorts';

export function dynamicComparedToBaseline(
    dataframe: DG.DataFrame,
    resCol: string,
    baselineVisitNum: string,
    blVisitColumn: string,
    newColName: string) {
    const dfName = dataframe.name;
    let grouped = dataframe.groupBy([PERSON_ID, resCol])
      .where(`${blVisitColumn} = ${baselineVisitNum}`)
      .aggregate();
    grouped.getCol(resCol).name = `BL_${resCol}`;
    dataframe.join(grouped,
      [PERSON_ID], [PERSON_ID],
      dataframe.columns.names(), [`BL_${resCol}`], DG.JOIN_TYPE.LEFT, true);
    dataframe.columns.addNewFloat(newColName)
      .init((i) => {
        if (dataframe.getCol(resCol).isNone(i) || dataframe.getCol(`BL_${resCol}`).isNone(i)) {
          return null;
        } else {
          return (dataframe.get(resCol, i) - dataframe.get(`BL_${resCol}`, i)) / dataframe.get(`BL_${resCol}`, i);
        }
      });
    dataframe.columns.remove(`BL_${resCol}`);
    dataframe.name = dfName;
  }

  export function joinCohorts(
    dataframe: DG.DataFrame) {
    const dfName = dataframe.name;
    const extractFromCohorts = cohorts.cohortsPivoted.columns.names().filter(it => it !== PERSON_ID);
     dataframe.join(cohorts.cohortsPivoted,
      [PERSON_ID], [PERSON_ID],
      dataframe.columns.names(), extractFromCohorts, DG.JOIN_TYPE.LEFT, true); 
    dataframe.name = dfName;
  }

  export function convertColToInt(dataframe: DG.DataFrame, colName: string) {
        let cohortCol = dataframe.col(colName);
        dataframe.columns.addNewInt(`${colName}_1`).init(i => parseInt(cohortCol.get(i)));
        dataframe.columns.remove(colName);
        dataframe.col(`${colName}_1`).name = colName;
  }

  export function convertColToString(dataframe: DG.DataFrame, colName: string) {
    let cohortCol = dataframe.col(colName);
    dataframe.columns.addNewString(`${colName}_1`).init(i => cohortCol.get(i).toString());
    dataframe.columns.remove(colName);
    dataframe.col(`${colName}_1`).name = colName;
}