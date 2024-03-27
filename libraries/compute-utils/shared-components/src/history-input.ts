/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, EMPTY, Observable, Subject, from} from 'rxjs';
import {catchError, finalize, shareReplay, switchMap, tap} from 'rxjs/operators';
import $ from 'cash-dom';
import {historyUtils} from '../../history-utils';
import '../css/history-input.css';

class DatabaseService {
  static getHistoryRuns(funcName: string, includeParams = true, skipDfLoad = false): Observable<DG.FuncCall[]> {
    return from((async () => {
      const res = await historyUtils.pullRunsByName(
        funcName, [
          {author: grok.shell.user},
          // EXPLAIN WHY FUNC.PARAMS
        ],
        {order: 'started'},
        [...(includeParams ? ['func.params']: []), 'session.user', 'options'],
        skipDfLoad,
      );
      return res;
    })());
  }
}

export const PICK_COLUMN_NAME = 'Pick' as const;
export const ID_COLUMN_NAME = 'ID' as const;

export abstract class HistoryInputBase<T = DG.FuncCall> extends DG.InputBase<T | null> {
  /**
   * Emitted when parictular historical run is chosen
   * @deprecated. Use {@link onInput} or {@link onChanged} instead
   * */
  public onHistoricalRunChosen = new BehaviorSubject<DG.FuncCall | null>(null);

  // Toggle to force historical runs udpate
  public experimentRunsUpdate = new Subject<boolean>();

  private store = {
    experimentRuns: new Observable<DG.FuncCall[]>(),
    experimentRunsDf: new BehaviorSubject(DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Pick', []),
      ...Object.keys(this._visibleColumnsForGrid).map((col) => DG.Column.fromStrings(col, [])),
    ])),
    isExperimentRunsLoading: new Subject<boolean>(),
    experimentRunData: new BehaviorSubject<string>(''),
  };
  private _historyGrid = DG.Viewer.grid(this.store.experimentRunsDf.value, {showRowHeader: false, showColumnGridlines: false, allowEdit: false});
  private _historyFilters = DG.Viewer.filters(this.store.experimentRunsDf.value, {title: 'Filters'});
  private _historyDialog = this.getHistoryDialog();

  private _visibleInput = ui.stringInput(this.label, '', null);
  private _visibleIcon = ui.iconFA('search', () => this.showSelectionDialog());
  private _chosenRun: DG.FuncCall | null = null;
  private _chosenRunId: string | null = null;

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
    private _visibleColumnsForFilter: string[] = [],
    // Load input/outputs
    private includeParams = true,
    // Load input/output dataframes
    private skipDfLoad = false,
  ) {
    const primaryInput = ui.stringInput(label, '', null);
    super(primaryInput.dart);

    this._visibleInput = primaryInput;
    this._visibleInput.readOnly = true;
    this._visibleInput.input.addEventListener('click', () => this.showSelectionDialog());
    this._visibleInput.addOptions(this._visibleIcon);

    this.store.experimentRuns = this.experimentRunsUpdate.pipe(
      tap(() => this.toggleLoaderExpRuns(true)),
      switchMap(() => DatabaseService.getHistoryRuns(this._funcName, this.includeParams, this.skipDfLoad).pipe(
        catchError((e) => {
          console.error(e);
          return EMPTY;
        }),
        finalize(() => this.toggleLoaderExpRuns(false)),
      )),
      shareReplay(1),
    );

    this.store.experimentRuns.subscribe((newRuns) => {
      const newRunsGridDf = DG.DataFrame.fromColumns([
        DG.Column.bool(PICK_COLUMN_NAME, newRuns.length).init(false),
        ...Object.entries(this._visibleColumnsForGrid).map((entry) => DG.Column.fromStrings(entry[0], newRuns.map(entry[1]))),
        DG.Column.fromStrings(ID_COLUMN_NAME, newRuns.map((newRun) => newRun.id)),
      ]);
      newRunsGridDf.onCurrentRowChanged.subscribe(() => {
        newRunsGridDf.getCol(PICK_COLUMN_NAME).init(false);
        newRunsGridDf.set(PICK_COLUMN_NAME, newRunsGridDf.currentRow.idx, true);
        const newId = newRunsGridDf.getCol(ID_COLUMN_NAME).get(newRunsGridDf.currentRow.idx);

        (this._historyDialog.getButton('OK') as HTMLButtonElement).disabled = false;
        this.setValue(newRuns.find((run) => run.id === newId) ?? null);
      });

      this.store.experimentRunsDf.next(newRunsGridDf);
    });

    this.store.experimentRunsDf.subscribe(() => {
      this.renderGridAndFilters();
    });

    this.store.isExperimentRunsLoading.subscribe((newValue) => {
      ui.setUpdateIndicator(this._visibleInput.root, newValue);
      ui.setUpdateIndicator(this._historyGrid.root, newValue);
    });

    this.experimentRunsUpdate.next();
  }

  private setCurrentRow() {
    if (this.getValue() || this.chosenRunId) {
      for (const row of this.store.experimentRunsDf.value.rows) {
        if (row.get(ID_COLUMN_NAME) === (this.getValue()?.id ?? this.chosenRunId)) {
          this.store.experimentRunsDf.value.currentRowIdx = row.idx;
          break;
        }
      }
    }
  }

  public showSelectionDialog() {
    // Bug: should re-render all viewers inside of the dialog
    this.renderGridAndFilters();

    (this._historyDialog.getButton('OK') as HTMLButtonElement).disabled = (this.store.experimentRunsDf.value.currentRow.idx === -1);

    this._historyDialog.show({
      modal: true,
      fullScreen: true,
      center: true,
      width: Math.max(400, window.screen.width - 200),
      height: Math.max(200, window.screen.height - 200),
    });
  }

  public renderGridAndFilters() {
    const newHistoryFilters = DG.Viewer.filters(this.store.experimentRunsDf.value, {title: 'Filters'});
    this._historyFilters.root.replaceWith(newHistoryFilters.root);
    this._historyFilters = newHistoryFilters;

    const newHistoryGrid = DG.Viewer.grid(this.store.experimentRunsDf.value, {showRowHeader: false, showColumnGridlines: false, allowEdit: false});
    this._historyGrid.root.replaceWith(newHistoryGrid.root);
    this._historyGrid = newHistoryGrid;

    this.styleHistoryGrid();
    this.styleHistoryFilters();

    this.setCurrentRow();
  }

  private getHistoryDialog() {
    const historyDialog = ui.dialog();

    historyDialog.onOK(async () => {
      if (this.getValue()) {
        let funcCall = this.getValue();
        if (!this.includeParams)
          funcCall = await historyUtils.loadRun(funcCall!.id);
        this.setValue(funcCall);
      }
      this.fireInput();
    });
    historyDialog.addButton('Refresh', () => {
      this.experimentRunsUpdate.next();
    }, 2, 'Reload list of runs');
    historyDialog.addButton('Discard', () => {
      this.setValue(null);
      this.historyGrid.dataFrame.currentRowIdx = -1;
      this.fireInput();
      historyDialog.close();
    }, 1, 'Discard the previous choice');
    $(historyDialog.getButton('CANCEL')).hide();

    historyDialog.add(
      ui.divH([
        ui.block([this._historyGrid.root]),
        ui.block([this._historyFilters.root], {style: {'margin-left': '10px', 'width': '300px', 'height': '100%'}}),
      ], {style: {'height': '100%'}}),
    );

    return historyDialog;
  }

  private styleHistoryFilters() {
    this._historyFilters.setOptions({columnNames: this._visibleColumnsForFilter});
    this._historyFilters.root.style.height = '100%';
    $(this._historyFilters.root).find('.d4-filter-group-header').hide();

    if (!this._visibleColumnsForFilter.length) $(this._historyFilters.root.parentElement).hide();
  }

  private styleHistoryGrid() {
    this._historyGrid.columns.byName(PICK_COLUMN_NAME)!.cellType = 'html';
    this._historyGrid.columns.byName(PICK_COLUMN_NAME)!.width = 30;
    this._historyGrid.root.style.height = '100%';
    this._historyGrid.setOptions({
      'showCurrentCellOutline': false,
      'allowEdit': false,
      'allowBlockSelection': false,
    });
    for (let i = 0; i < this._historyGrid.columns.length; i++) {
      const col = this._historyGrid.columns.byIndex(i)?.column;
      if (col && col.type === DG.TYPE.DATE_TIME)
        this._historyGrid.columns.byIndex(i)!.format = 'MMM d yyyy HH:mm';
    }

    this._historyGrid.onCellPrepare((cell) => {
      if (cell.tableColumn?.name === PICK_COLUMN_NAME) {
        if (cell.isColHeader)
          cell.customText = '';

        cell.element = cell.cell.value ? ui.div([ui.iconFA('check', null, 'Selected')], {style: {'text-align': 'center', 'margin': '6px'}}): ui.div();
      }
    });

    this._historyGrid.columns.setVisible([PICK_COLUMN_NAME, ...Object.keys(this._visibleColumnsForGrid).map((colName) => colName)]);
  }

  protected setValue(val: DG.FuncCall | null) {
    this._chosenRun = val;
    this._chosenRunId = this._chosenRun?.id ?? null;
    this.onHistoricalRunChosen.next(val);
    this._visibleInput.value = val ? this._stringValueFunc(val): 'No run chosen';
    if (!this.notify) return;

    this.fireChanged();
  }

  protected getValue() {
    return this._chosenRun;
  }

  get chosenRun() {
    return this._chosenRun;
  }

  get chosenRunId() {
    return this._chosenRunId;
  }

  set chosenRunId(id: string | null) {
    this._chosenRunId = id;
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

  set stringValue(val: string) {
    this._visibleInput.value = val;
  }

  get stringValue() {
    return this._visibleInput.value;
  }

  get iconRoot() {
    return this._visibleIcon;
  }

  public toggleLoaderExpRuns(newState: boolean) {
    this.store.isExperimentRunsLoading.next(newState);
  }
}

export class HistoryInput extends HistoryInputBase {
  set value(val: DG.FuncCall | null) {
    this.setValue(val);
  }

  get value() {
    return this.getValue();
  }
}

export class HistoryInputJSON extends HistoryInputBase<string | null> {
  get value() {
    const funccall = this.getValue();
    if (funccall) {
      return JSON.stringify({
        visualValue: this.stringValue,
        // BUG HERE: if GUID is present in the string, DB fails to save it
        id: funccall.id.split('-'),
      });
    } else
      return null;
  }

  set value(jsonVal: string | null) {
    this.setValue(null);
    this.stringValue = jsonVal ? JSON.parse(jsonVal).visualValue: 'No run is chosen';
    // Recovering initial GUID here
    this.chosenRunId = jsonVal ? JSON.parse(jsonVal).id.join('-'): '';
  }
}
