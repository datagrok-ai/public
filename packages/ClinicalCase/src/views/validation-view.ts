import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { validationRulesList } from "../package";
import { pinnacleRuleIdColumnName, validationResultRuleIdColumn } from "../validation/constants";
import { createValidationDataFrame } from '../validation/validation-utils';

export class ValidationView extends DG.ViewBase {

  resultsDataframe: DG.DataFrame;
  rulesDataframe: DG.DataFrame;
  resultsGrid: DG.Grid;
  rulesGrid: DG.Grid;

  constructor(domain?: string) {
    super();

    this.resultsDataframe = domain? this.createFilteredDataframe(study.validationResults, domain) :study.validationResults;
    this.rulesDataframe = validationRulesList;

    this.resultsGrid = this.resultsDataframe.plot.grid();
    this.rulesGrid = this.rulesDataframe.plot.grid();

    grok.data.linkTables(this.rulesDataframe, this.resultsDataframe,
      [ `${pinnacleRuleIdColumnName}` ], [ `${validationResultRuleIdColumn}` ],
      [ DG.SYNC_TYPE.CURRENT_ROW_TO_FILTER, DG.SYNC_TYPE.SELECTION_TO_SELECTION ]);


    let uniqueViolatedRuleIds = this.getUniqueErrorIds(this.resultsDataframe);

    this.rulesGrid.onCellPrepare(function (gc) {
      if (!gc.isTableCell)
        return;
      if (gc.gridColumn.idx === 1 && uniqueViolatedRuleIds.has(gc.cell.value))
        gc.style.backColor = 0xFFFF0000;
    });

    this.generateUI(this.resultsGrid, this.rulesGrid);
  }

  private generateUI(resultsGrid: any, rulesGrid: any) {
    this.root.appendChild(ui.div([
      ui.divH(
        [
          ui.block([ ui.h1('List of rules'), rulesGrid ], { style: { width: '100%' } }),
        ]),
      ui.divH([
        ui.block([ ui.h1('Detected errors'), resultsGrid ]),
      ])
    ]));
  }

  private getUniqueErrorIds(resultsDataframe: DG.DataFrame) {
    const uniqueIds = new Set();
    let column = resultsDataframe.columns.byName(validationResultRuleIdColumn);
    let rowCount = resultsDataframe.rowCount;
    for (let i = 0; i < rowCount; i++)
      uniqueIds.add(column.get(i));
    return uniqueIds;
  }

  private createFilteredDataframe(resultsDataframe: DG.DataFrame, domain: string){
    const res = createValidationDataFrame();
    let column = resultsDataframe.columns.byName('Domain');
    let rowCount = resultsDataframe.rowCount;
    for (let i = 0; i < rowCount; i++){
      if(column.get(i) === domain){
        const array = [];
        for (let col of resultsDataframe.columns)
            array.push(col.get(i));
        res.rows.addNew(array);
      }
    }
    return res;
  }

}

