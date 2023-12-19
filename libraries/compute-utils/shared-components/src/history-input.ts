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
  static getHistoryRuns(funcName: string, includeParams = true): Observable<DG.FuncCall[]> {
    return from((async () => {
      const res = await historyUtils.pullRunsByName(
        funcName, [
          {author: grok.shell.user},
          // EXPLAIN WHY FUNC.PARAMS
        ], {order: 'started'}, [...(includeParams ? ['func.params']: []), 'session.user', 'options'],
      );
      return res;
    })());
  }
}

const PICK_COLUMN_NAME = 'Pick' as const;
const ID_COLUMN_NAME = 'ID' as const;

export class HistoryInput extends DG.InputBase<DG.FuncCall | null> {
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
  private _value: DG.FuncCall | null = null;

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
    // Load input/output dataframes
    private includeParams = true,
  ) {
    const primaryInput = ui.stringInput(label, '', null);
    super(primaryInput.dart);

    this._visibleInput = primaryInput;
    this._visibleInput.readOnly = true;
    this._visibleInput.input.addEventListener('click', () => this.showSelectionDialog());
    this._visibleInput.addOptions(ui.iconFA('search', () => this.showSelectionDialog()));

    this.store.experimentRuns = this.experimentRunsUpdate.pipe(
      tap(() => this.toggleLoaderExpRuns(true)),
      switchMap(() => DatabaseService.getHistoryRuns(this._funcName, this.includeParams).pipe(
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
        this._value = newRuns.find((run) => run.id === newId) ?? null;
      });
      this.store.experimentRunsDf.next(newRunsGridDf);
    });

    this.store.experimentRunsDf.subscribe((newDf) => {
      // Bug: should re-render all viewers inside of the dialog
      const newHistoryFilters = DG.Viewer.filters(newDf, {title: 'Filters'});
      this._historyFilters.root.replaceWith(newHistoryFilters.root);
      this._historyFilters = newHistoryFilters;

      const newHistoryGrid = DG.Viewer.grid(newDf, {showRowHeader: false, showColumnGridlines: false, allowEdit: false});
      this._historyGrid.root.replaceWith(newHistoryGrid.root);
      this._historyGrid = newHistoryGrid;

      this.styleHistoryGrid();
      this.styleHistoryFilters();
    });

    this.store.isExperimentRunsLoading.subscribe((newValue) => {
      ui.setUpdateIndicator(this._visibleInput.root, newValue);
    });

    this.experimentRunsUpdate.next();
  }

  public showSelectionDialog() {
    this._historyDialog.show({
      modal: true,
      fullScreen: true,
      center: true,
      width: Math.max(400, window.screen.width - 200),
      height: Math.max(200, window.screen.height - 200),
    });
  }

  private getHistoryDialog() {
    const historyDialog = ui.dialog();

    historyDialog.onOK(async () => {
      if (this.value) {
        let funcCall = this.value;
        if (!this.includeParams)
          funcCall = await historyUtils.loadRun(funcCall.id);
        this.value = funcCall;
      }
      this.fireInput();
    });
    historyDialog.addButton('Refresh', () => {
      this.experimentRunsUpdate.next();
    }, 2, 'Reload list of runs');
    historyDialog.addButton('Discard', () => {
      this.value = null;
      this.historyGrid.dataFrame.currentRowIdx = -1;
      historyDialog.close();
    }, 1, 'Discard the previous choice');
    $(historyDialog.getButton('CANCEL')).hide();
    (historyDialog.getButton('OK') as HTMLButtonElement).disabled = (this.store.experimentRunsDf.value.currentRow.idx === -1);

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
    this._value = val;
    this.onHistoricalRunChosen.next(val);
    this._visibleInput.value = val ? this._stringValueFunc(val): '';
    if (!this.notify) return;

    this.fireChanged();
  }

  get value() {
    return this._value;
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
