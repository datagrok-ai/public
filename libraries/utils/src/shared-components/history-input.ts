/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, EMPTY, Observable, Subject, from} from 'rxjs';
import {catchError, finalize, shareReplay, switchMap, tap, takeUntil} from 'rxjs/operators';
import $ from 'cash-dom';
import {historyUtils} from '../history-utils';

class DatabaseService {
  static getHistoryRuns(funcName: string): Observable<DG.FuncCall[]> {
    return from((async () => {
      const res = await historyUtils.pullRunsByName(
        funcName, [
          {isShared: true}, {author: grok.shell.user},
        ], {order: 'started'}, ['func.params', 'session.user', 'options']
      );
      return res;
    })());
  }
}

export class HistoryInput {
  // emitted when parictular historical run is chosen
  public onHistoricalRunChosen = new BehaviorSubject<DG.FuncCall | null>(null);

  // Toggle to force historical runs udpate
  public experimentRunsUpdate = new Subject<boolean>();

  private store = {
    experimentRuns: new Observable<DG.FuncCall[]>(),
    experimentRunsDf: new BehaviorSubject(DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Pick', []),
      ...Object.keys(this._visibleColumnsForGrid).map((col) => DG.Column.fromStrings(col, []))
    ])),
    isExperimentRunsLoading: new Subject<boolean>(),
    experimentRunData: new BehaviorSubject<string>(''),
  };
  private _historyGrid = DG.Viewer.grid(this.store.experimentRunsDf.value, {showRowHeader: false, showColumnGridlines: false, allowEdit: false});
  private _historyFilters = DG.Viewer.filters(this.store.experimentRunsDf.value, {title: 'Filters'});

  private _visibleInput = ui.stringInput(this.label, '', null);

  constructor(
    // Label placed before next to the input
    private label: string,
    // Function to look history for
    private _funcName: string,
    // String value to be displayed in the input field
    private _stringValueFunc: (currentRun: DG.FuncCall) => string,
    // Map of grid columns and functions to retrieve their actual values from DG.FuncCall
    private _visibleColumnsForGrid: Record<string, (currentRun: DG.FuncCall) => string>,
    // Array of grid columns visible in the filter viever
    private _visibleColumnsForFilter: string[] = []
  ) {
    this._historyGrid.columns.byName('Pick')!.cellType = 'html';
    this._historyGrid.columns.byName('Pick')!.width = 30;

    this._visibleInput.input.addEventListener('click', () => {
      const _historyDialog = ui.dialog();

      $(_historyDialog.root).addClass('arrow-box');

      _historyDialog.onOK(() => { });
      _historyDialog.onCancel(() => { this.onHistoricalRunChosen.next(null); });

      $(_historyDialog.root).find('.d4-dialog-header').hide();
      _historyDialog.root.style.boxShadow = '';
      $(_historyDialog.root).removeClass('ui-form');

      const newHistoryFilters = DG.Viewer.filters(this.store.experimentRunsDf.value, {title: 'Filters'});
      this._historyFilters.root.replaceWith(newHistoryFilters.root);
      this._historyFilters = newHistoryFilters;

      const newHistoryGrid = DG.Viewer.grid(this.store.experimentRunsDf.value, {showRowHeader: false, showColumnGridlines: false, allowEdit: false});
      this._historyGrid.root.replaceWith(newHistoryGrid.root);
      this._historyGrid = newHistoryGrid;

      this._historyGrid.columns.byName('Pick')!.cellType = 'html';
      this._historyGrid.columns.byName('Pick')!.width = 30;

      this._historyGrid.onCellPrepare((cell) => {
        if (cell.tableColumn?.name === 'Pick') {
          if (cell.isColHeader)
            cell.customText = '';

          if (cell.cell.value === true)
            cell.element = ui.div([ui.iconFA('check')], {style: {'text-align': 'center', 'margin': '6px'}});
          else
            cell.element = ui.div();
        }
      });

      this._historyGrid.columns.setVisible(['Pick', ...Object.keys(this._visibleColumnsForGrid).map((colName) => colName)]);
      this._historyFilters.setOptions({columnNames: this._visibleColumnsForFilter});
      $(this._historyFilters.root).find('.d4-filter-group-header').hide();

      if (!this._visibleColumnsForFilter.length) $(this._historyFilters.root.parentElement).hide();

      _historyDialog.add(
        ui.divH([
          ui.block([this._historyGrid.root]),
          ui.block([this._historyFilters.root], {style: {'margin-left': '10px', 'width': '225px'}}),
        ])
      );

      _historyDialog.show({
        x: this._visibleInput.input.getBoundingClientRect().x + this._visibleInput.input.clientWidth + 10,
        y: this._visibleInput.input.getBoundingClientRect().y - 30,
        width: Math.max(this._historyGrid.root.clientWidth + this._historyFilters.root.clientWidth + 30, 600),
      });
    });
    this._visibleInput.readOnly = true;

    this.store.experimentRuns = this.experimentRunsUpdate.pipe(
      tap(() => this.toggleLoaderExpRuns(true)),
      switchMap(() => DatabaseService.getHistoryRuns(this._funcName).pipe(
        catchError(() => EMPTY),
        finalize(() => this.toggleLoaderExpRuns(false))
      )),
      shareReplay(1)
    );

    this.store.experimentRuns.pipe(takeUntil(from(new Subject<boolean>()))).subscribe((newRuns) => {
      const newRunsGridDf = DG.DataFrame.fromColumns([
        DG.Column.bool('Pick', newRuns.length).init(false),
        ...Object.entries(this._visibleColumnsForGrid).map((entry) => DG.Column.fromStrings(entry[0], newRuns.map(entry[1]))),
        DG.Column.fromStrings('ID', newRuns.map((newRun) => newRun.id))
      ]);
      newRunsGridDf.onCurrentRowChanged.subscribe(() => {
        newRunsGridDf.getCol('Pick').init(false);
        newRunsGridDf.set('Pick', newRunsGridDf.currentRow.idx, true);
        const newId = newRunsGridDf.getCol('ID').get(newRunsGridDf.currentRow.idx);
        const chosenRun = newRuns.find((run) => run.id === newId);
        if (chosenRun)
          this.onHistoricalRunChosen.next(chosenRun);
      });
      this.store.experimentRunsDf.next(newRunsGridDf);
    });

    this.store.experimentRunsDf.subscribe((newDf) => {
      this.historyGrid.dataFrame = newDf;
      this.historyFilters.dataFrame = newDf;
    });

    this.onHistoricalRunChosen.subscribe((val) => {
      if (val)
        this._visibleInput.value = this._stringValueFunc(val);
      else
        this._visibleInput.value = '';
    });

    this.store.isExperimentRunsLoading.subscribe((newValue) => {
      ui.setUpdateIndicator(this._visibleInput.root, newValue);
    });
  }

  // Viever object of grid
  get historyGrid() {
    return this._historyGrid;
  }

  // Viever object of filters
  get historyFilters() {
    return this._historyFilters;
  }

  // HTML root of component
  get root() {
    return this._visibleInput.root;
  }

  set value(val: DG.FuncCall | null) {
    if (val)
      this._visibleInput.value = this._stringValueFunc(val);
    else
      this._visibleInput.value = '';
  }

  get value() {
    return this.onHistoricalRunChosen.value;
  }

  set stringValue(val: string) {
    this._visibleInput.value = val;
  }

  get stringValue() {
    return this._visibleInput.value;
  }

  public toggleLoaderExpRuns(newState: boolean) {
    this.store.isExperimentRunsLoading.next(newState);
  }
}
