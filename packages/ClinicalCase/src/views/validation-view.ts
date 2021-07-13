import { ClinicalCaseView } from "../clinical-case-view";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { validationRulesList } from "../package";
import { pinnacleRuleIdColumnName } from "../validation/constants";

export class ValidationView extends DG.ViewBase {

  constructor() {
    super();

    let resultsGrid = study.validationResults.plot.grid();
    let rulesGrid = validationRulesList.plot.grid();

    validationRulesList.onCurrentRowChanged.subscribe(args => {
      const currentRowIdx = validationRulesList.currentRow.idx;
      const currentRuleId = validationRulesList.get(pinnacleRuleIdColumnName, currentRowIdx);
     /*  study.validationResults.rows.filter((row) => row.get('Violated rule id') === currentRuleId);
      resultsGrid.invalidate(); */
    });

    this.root.appendChild(ui.div([
      ui.divH([ rulesGrid ], { style: { width: '100%' } }),
      ui.divH([ resultsGrid ], { style: { width: '100%' } })
    ]));
    
  }
}