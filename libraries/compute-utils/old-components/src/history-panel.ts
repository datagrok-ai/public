/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import {historyUtils} from '../../history-utils';
import '../css/history-panel.css';
import {HistoricalRunsList} from './history-list';
import {filter} from 'rxjs/operators';

export class HistoryPanel extends DG.Widget {
  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCalls are called for comparison. Contains FuncCalls' IDs
  public onComparison = new Subject<string[]>();

  // Emitted when FuncCall is edited
  public afterRunEdited = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is deleted
  public afterRunDeleted = new Subject<DG.FuncCall>();

  public allRunsFetch = new Subject<true>();

  private historyList = new HistoricalRunsList([], {fallbackText: 'No runs are found in history',
    showActions: true, showBatchActions: true, isHistory: true});
  private panel = this.historyList.root;

  private updateHistoryPane(historyRuns: DG.FuncCall[]) {
    this.historyList.updateRuns(historyRuns);
  };

  public showEditDialog(editedCall: DG.FuncCall) {
    this.historyList.editRun(editedCall);
  }

  public addRun(newRun: DG.FuncCall) {
    this.historyList.addRun(newRun);
  }

  public get history() {
    return this.historyList.runs;
  }

  public get onCompactModeChanged() {
    return this.historyList.onCompactModeChanged;
  }

  public get compactMode() {
    return this.historyList.compactMode;
  }

  constructor(
    private func: DG.Func,
  ) {
    super(ui.box(ui.divText('No historical runs loaded', 'hp-no-elements-label'), {style: {'width': '100%', 'height': '100%'}}));

    const clickedSub = this.historyList.onChosen
      .pipe(
        filter((call) => !!call),
      ).subscribe((clickedCall) => this.onRunChosen.next(clickedCall!.id));

    const comparisonSub = this.historyList.onComparisonCalled.subscribe((ids) => this.onComparison.next(ids));

    const allRunsFetch = this.allRunsFetch.subscribe(async () => {
      ui.setUpdateIndicator(this.root, true);
      historyUtils.pullRunsByName(this.func.name, [{author: grok.shell.user}], {}, ['session.user', 'options'])
        .then((historicalRuns) => {
          ui.empty(this.root);
          this.root.appendChild(this.panel);
          this.updateHistoryPane(historicalRuns);
        })
        .catch((e) => grok.shell.error(e))
        .finally(() => ui.setUpdateIndicator(this.root, false));
    });

    const onMetadataEdit = this.historyList.onMetadataEdit.subscribe((editedCall) => {
      this.afterRunEdited.next(editedCall);
    });

    const onDeleteSub = this.historyList.onDelete.subscribe((deleteCall) => {
      this.afterRunDeleted.next(deleteCall);
    });

    this.subs.push(
      clickedSub,
      allRunsFetch,
      onMetadataEdit,
      comparisonSub,
      onDeleteSub,
    );

    this.allRunsFetch.next();
  }
}
