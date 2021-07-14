import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { validationRulesList } from "../package";
import { pinnacleRuleIdColumnName, validationResultRuleIdColumn } from "../validation/constants";

export class ValidationView extends DG.ViewBase {

  constructor() {
    super();

    let resultsGrid = study.validationResults.plot.grid();
    let rulesGrid = validationRulesList.plot.grid();

    let uniqueViolatedRuleIds = this.getUniqueErrorIds();

    validationRulesList.onCurrentRowChanged.subscribe(args => {
      const currentRowIdx = validationRulesList.currentRow.idx;
      const currentRuleId = validationRulesList.get(pinnacleRuleIdColumnName, currentRowIdx);
      study.validationResults.rows.match({ 'Violated rule id': `${currentRuleId}`}).filter();
      resultsGrid.invalidate();
    });

    rulesGrid.onCellPrepare(function (gc) {
      if (!gc.isTableCell)
        return;
      if (gc.gridColumn.idx === 1 && uniqueViolatedRuleIds.has(gc.cell.value))
        gc.style.backColor = 0xFFFF0000;    
    });

    this.generateUI(resultsGrid, rulesGrid);
  }

  private generateUI(resultsGrid: any, rulesGrid: any){
    this.root.appendChild(ui.div([
      ui.divH(
        [
          ui.block([ui.h2('List of rules'), rulesGrid], { style: { width: '100%' } }),
        ]),
      ui.divH([
        ui.block([ui.h2('Detected errors'), resultsGrid]),
      ])
    ]));
  }

  private getUniqueErrorIds(){
    const uniqueIds = new Set();
    let column = study.validationResults.columns.byName(validationResultRuleIdColumn);
    let rowCount = study.validationResults.rowCount;
    for (let i = 0; i < rowCount; i++)
      uniqueIds.add(column.get(i));
    return uniqueIds;
  }

}

