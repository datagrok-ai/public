import * as grok from 'datagrok-api/grok';
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import { study } from "../clinical-study";
import { validationRulesList } from "../package";
import { pinnacleRuleIdColumnName, validationResultRuleIdColumn } from "../validation/constants";
import { createRulesDataFrame, createValidationDataFrame } from '../validation/validation-utils';
import { View } from 'datagrok-api/dg';

export class ValidationView extends DG.ViewBase {

  resultsDataframe: DG.DataFrame;
  rulesDataframe: DG.DataFrame;
  resultsGrid: DG.Grid;
  rulesGrid: DG.Grid;
  errorsByDomain: any;
  domains: any;

  constructor(errorsMap: any) {
    super();

    this.errorsByDomain = errorsMap;
    this.resultsDataframe = study.validationResults;
    this.domains = study.domains;

    let uniqueViolatedRuleIds = this.getUniqueErrorIds(study.validationResults);
    this.rulesDataframe = this.getViolatedRulesDataframe(validationRulesList, uniqueViolatedRuleIds);

    //this.resultsGrid = this.resultsDataframe.plot.grid();
    this.rulesGrid = this.rulesDataframe.plot.grid();

    grok.data.linkTables(this.rulesDataframe, this.resultsDataframe,
      [ `${pinnacleRuleIdColumnName}` ], [ `${validationResultRuleIdColumn}` ],
      [ DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION, DG.SYNC_TYPE.CURRENT_ROW_TO_SELECTION ]);



    // this.rulesGrid.onCellPrepare(function (gc) {
    //   if (!gc.isTableCell)
    //     return;
    //   if (gc.gridColumn.idx === 1 && uniqueViolatedRuleIds.has(gc.cell.value))
    //     gc.style.backColor = 0xFFFF0000;
    // });

    this.generateUI(this.resultsGrid, this.rulesGrid);
  }

  private generateUI(resultsGrid: any, rulesGrid: any) {

    this.root.appendChild(
      ui.splitV([

        ui.box(ui.panel([
          ui.h1('Violated Rules')
        ]), { style: { maxHeight: '30px' } }),

        DG.Viewer.grid(this.rulesDataframe).root,

        ui.box(ui.panel([
          ui.h1('Errors')
        ]), { style: { maxHeight: '30px' } }),

        ui.tabControl(this.createErrorsTabControl()).root
      ], { style: { width: '100%', height: '100%' } })
    );

  }

  private createErrorsTabControl() {
    const tabControl = {}
    Object.keys(this.errorsByDomain).forEach(key => {
      const domainDataframe = this.domains[ key ];
      this.addRowNumberColumn(domainDataframe, key);
      tabControl[`${key.toUpperCase()} (${this.errorsByDomain[ key ]})`] = DG.Viewer.grid(domainDataframe).root;
      grok.data.linkTables(this.resultsDataframe, domainDataframe,
        [ `Row number`,  'Domain'], [ `Row number`,  'Domain lower case' ],
        [ DG.SYNC_TYPE.SELECTION_TO_SELECTION, DG.SYNC_TYPE.SELECTION_TO_FILTER ]);
    })
    return tabControl;
  }


  private addRowNumberColumn(dataframe: DG.DataFrame, domain: string){
    if(!dataframe.columns.contains('Row number')){
      dataframe.columns.addNewInt('Row number').init((i) => i+1);
    }
    if(!dataframe.columns.contains('Domain lower case')){
      dataframe.columns.addNewString('Domain lower case').init((i) => domain);
    }
  }

  private getUniqueErrorIds(resultsDataframe: DG.DataFrame) {
    const uniqueIds = new Set();
    let column = resultsDataframe.columns.byName(validationResultRuleIdColumn);
    let rowCount = resultsDataframe.rowCount;
    for (let i = 0; i < rowCount; i++)
      uniqueIds.add(column.get(i));
    return uniqueIds;
  }

  private getViolatedRulesDataframe(rules: DG.DataFrame, uniqueViolatedRuleIds: any) {
    const res = createRulesDataFrame();
    let column = rules.columns.byName(pinnacleRuleIdColumnName);
    let rowCount = rules.rowCount;
    for (let i = 0; i < rowCount; i++) {
      if (uniqueViolatedRuleIds.has(column.get(i))) {
        const array = [];
        for (let col of rules.columns)
          array.push(col.get(i));
        res.rows.addNew(array);
      }
    }
    return res;
  }


  private createFilteredDataframe(resultsDataframe: DG.DataFrame, domain: string) {
    const res = createValidationDataFrame();
    let column = resultsDataframe.columns.byName('Domain');
    let rowCount = resultsDataframe.rowCount;
    for (let i = 0; i < rowCount; i++) {
      if (column.get(i) === domain) {
        const array = [];
        for (let col of resultsDataframe.columns)
          array.push(col.get(i));
        res.rows.addNew(array);
      }
    }
    return res;
  }

}

