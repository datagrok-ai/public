import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export class SubstituentAnalysisViewer extends DG.JsViewer {
  initialized: boolean = false;
  name = 'substituent';
  grouppingCols: string[] = [];
  colsToAnalyze: string[] = [];
  totalCols: string[] = [];
  grouppedDf: DG.DataFrame | undefined = undefined;
  grouppedGridDiv = ui.div();
  mainView = ui.splitV([]);
  historgamMode: DG.InputBase<boolean> | undefined = undefined;

  constructor() {
    super();
  }

  init(): void {
    this.initialized = true;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  async onTableAttached(): Promise<void> {
    this.init();
    const grouppingColsDiv = ui.div();
    const colsToAnalyzeDiv = ui.div();
    this.totalCols = this.dataFrame.columns.names().concat(['']);
    this.updateColumnChoices(this.grouppingCols, 'Group by', grouppingColsDiv);
    this.updateColumnChoices(this.colsToAnalyze, 'Values', colsToAnalyzeDiv);
    this.mainView.append(ui.box(ui.panel([grouppingColsDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(ui.box(ui.panel([colsToAnalyzeDiv], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.historgamMode = ui.switchInput('Historgam', true, () => {
      this.updateGrid();
    });
    this.mainView.append(
      ui.box(ui.panel([this.historgamMode.root], {style: {padding: '0px'}}), {style: {maxHeight: '30px'}}));
    this.mainView.append(this.grouppedGridDiv);
    this.root.append(this.mainView);
  }

  updateColumnChoices(selectedCols: string[], choicesInputsName: string, choicesDiv: HTMLDivElement) {
    const colsChoicesDiv = ui.divH([]);
    colsChoicesDiv.append(ui.divText(
      choicesInputsName, {style: {fontWeight: 'bold', paddingTop: '10px', paddingRight: '10px'}}));
    for (let i = 0; i < selectedCols.length + 1; i++) {
      const selectedValue = selectedCols.length === i ? '' : selectedCols[i];
      const groupChoiceInput = ui.choiceInput('', selectedValue, this.totalCols);
      ui.tooltip.bind(groupChoiceInput.root, () => selectedValue);
      groupChoiceInput.onChanged(() => {
        if (groupChoiceInput.value === '')
          selectedCols.splice(i, 1);
        else {
          selectedCols.length === i ? selectedCols.push(groupChoiceInput.value!) :
            selectedCols[i] = groupChoiceInput.value!;
        }
        this.updateColumnChoices(selectedCols, choicesInputsName, choicesDiv);
        this.updateGrid();
      });
      groupChoiceInput.input.style.width = '100px';
      colsChoicesDiv.append(groupChoiceInput.root);
    }
    ui.empty(choicesDiv);
    choicesDiv.append(colsChoicesDiv);
  }

  updateGrid() {
    this.grouppedDf = this.dataFrame.groupBy(this.grouppingCols).aggregate();
    const grid = this.grouppedDf.plot.grid();
    grid.root.style.width = '100%';
    grid.root.style.height = '100%';
    for (const col of this.grouppingCols)
      grid.col(col)!.width = 150;

    for (const col of this.colsToAnalyze) {
      this.grouppedDf.columns.addNewString(col).init((i) => '');
      grid.col(col)!.cellType = 'html';
      grid.col(col)!.width = 150;
    }
    grid.onCellPrepare((gc) => {
      if (gc.isTableCell && gc.cell.value != null && this.colsToAnalyze.includes(gc.gridColumn.name)) {
        const condition: {[key: string]: any} = {};
        for (const col of this.grouppingCols)
          condition[col] = this.grouppedDf!.get(col, gc.tableRowIndex!);
        const df = this.dataFrame.clone(null, this.grouppingCols.concat([gc.gridColumn.name]));
        if (this.historgamMode!.value) {
          gc.style.element = df.plot
            .histogram({value: gc.gridColumn.name,
              showColumnSelector: false, showFilteredOutRows: true, filteringEnabled: false}).root;
          df.rows.match(condition).filter();
          console.log(`filter count: ${df.filter.trueCount}`);
          const test = df.rows.match(condition).toDataFrame();
          console.log(`df rowCount: ${test.rowCount}`);
        } else {
          df.columns.addNewString('group').init((j) => {
            for (const group of this.grouppingCols) {
              if (df.get(group, j) !== this.grouppedDf!.get(group, gc.tableRowIndex!))
                return 'other';
            }
            return 'current_group';
          });
          gc.style.element = gc.style.element = df.plot
            .histogram({value: gc.gridColumn.name, split: 'group', bins: 5,
              showColumnSelector: false, showFilteredOutRows: true, filteringEnabled: false}).root;
        }
      }
    });
    ui.empty(this.grouppedGridDiv);
    this.grouppedGridDiv.append(grid.root);
  }
}
