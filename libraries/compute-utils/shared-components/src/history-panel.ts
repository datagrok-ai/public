/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {BehaviorSubject, Subject, merge} from 'rxjs';
import {historyUtils} from '../../history-utils';
import '../css/history-panel.css';
import {getObservable, properUpdateIndicator} from '../../function-views/src/shared/utils';
import {HistoricalRunsList} from './history-list';

type FilterOptions = {
  isOnlyFavorites: boolean,
  text: string,
  author?: DG.User,
  startedAfter?: dayjs.Dayjs,
  tags?: string[]
}

class HistoryFilter extends DG.Widget {
  public onFilteringChanged = new BehaviorSubject<FilterOptions>({text: '', isOnlyFavorites: false});

  private onFuncCallListChanged = new BehaviorSubject<DG.FuncCall[]>([]);

  constructor() {
    super(ui.div());

    const favoritesInput = ui.switchInput('Only favorites', false);

    const textInput = ui.stringInput('Search', '', null, {clearIcon: true, placeholder: 'Type here to filter...'});
    const dateInput = ui.dateInput('Started after', dayjs().subtract(1, 'week'));
    dateInput.addPatternMenu('datetime');
    const tagInput = ui.choiceInput<string>('Tag', 'Choose tag to filter', ['Choose tag to filter']);

    const filterTagEditor = DG.TagEditor.create();
    const dummyInput = ui.stringInput(' ', '');
    dummyInput.input.replaceWith(filterTagEditor.root);
    ui.setDisplay(dummyInput.root, false);

    const addTagSub = getObservable<string | null>(tagInput.onChanged.bind(tagInput))
      .subscribe(() => {
        const tag = tagInput.value;
        //@ts-ignore
        if (!!tag && tag !== 'Choose tag to filter' && !filterTagEditor.tags.includes(tag))
          filterTagEditor.addTag(tag);
        tagInput.notify = false;
        tagInput.value = 'Choose tag to filter';
        tagInput.notify = true;
      });

    const filterChangeSub = merge(
      getObservable<boolean>(favoritesInput.onInput.bind(favoritesInput)),
      DG.debounce(getObservable<string>(textInput.onInput.bind(textInput))),
      getObservable<dayjs.Dayjs>(dateInput.onInput.bind(dateInput)),
      getObservable<string[]>(filterTagEditor.onChanged.bind(filterTagEditor)),
    ).subscribe(() => {
      const isOnlyFavorites = favoritesInput.value;
      const text = textInput.value;
      const startedAfter = dateInput.value ?? undefined;
      //@ts-ignore
      const tags: string[] = filterTagEditor.tags.filter((tag) => !!tag) ?? [];
      this.onFilteringChanged.next({isOnlyFavorites, text, startedAfter, tags});

      if (tags.length === 0)
        ui.setDisplay(dummyInput.root, false);
      else
        ui.setDisplay(dummyInput.root, true);
    });

    const updateInputsSub = this.onFuncCallListChanged.subscribe((newFuncCalls) => {
      properUpdateIndicator(tagInput.root, true);
      const newTags = newFuncCalls
        .map((run) => run.options['tags'] as string[])
        .reduce((acc, runTags) => {
          if (!!runTags) runTags.forEach((runTag) => {if (!acc.includes(runTag)) acc.push(runTag);});
          return acc;
        }, [] as string[]);

      tagInput.items = ['Choose tag to filter', ...newTags];
      dateInput.value = newFuncCalls.reduce((minDate, call) => {
        if (call.started < minDate)
          minDate = call.started;

        return minDate;
      }, dayjs()).subtract(1, 'week');
      properUpdateIndicator(tagInput.root, false);
    });

    this.subs.push(filterChangeSub, addTagSub, updateInputsSub);

    const form = ui.divV([
      favoritesInput,
      textInput,
      dateInput,
      tagInput,
      dummyInput,
    ], 'ui-form ui-form-wide ui-form-left');

    this.root.appendChild(form);
  }

  public addTag(funccall: DG.FuncCall) {
    const newFunccalls = this.onFuncCallListChanged.value.filter((fc) => fc.id !== funccall.id);
    this.onFuncCallListChanged.next([...newFunccalls, funccall]);
  }

  public removeTag(funccall: DG.FuncCall) {
    const newFunccalls = this.onFuncCallListChanged.value.filter((fc) => fc.id !== funccall.id);
    this.onFuncCallListChanged.next(newFunccalls);
  }

  public updateTagList(funccalls: DG.FuncCall[]) {
    this.onFuncCallListChanged.next(funccalls);
  }
}

export class HistoryPanel extends DG.Widget {
  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCalls are called for comparison. Contains FuncCalls' IDs
  public onComparison = new Subject<string[]>();

  // Emitted when FuncCall is edited
  public onRunEdited = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is deleted
  public afterRunDeleted = new Subject<DG.FuncCall>();

  public allRunsFetch = new Subject<true>();

  private historyFilter = new HistoryFilter();
  private historyList = new HistoricalRunsList([], {fallbackText: 'No runs are found in history',
    showAuthorIcon: true, showDelete: true, showEdit: true, showFavorite: true, showCompare: true});
  private panel = ui.divV([
    this.historyFilter.root,
    this.historyList.root,
  ]);
  private allRuns = new BehaviorSubject<DG.FuncCall[]>([]);
  private get historyRuns() {
    return this.allRuns.value.filter((run) => run.author.id === grok.shell.user.id);
  }
  private get filteredHistoryRuns() {
    const filteringOptions = this.historyFilter.onFilteringChanged.value;
    return this.historyRuns.filter((run) => {
      const isFavorite = filteringOptions.isOnlyFavorites ? run.options['isFavorite'] : true;
      const startedAfter = (filteringOptions.startedAfter ? run.started > filteringOptions.startedAfter : true);
      const titleContainsText = (!filteringOptions.text.length) ? true : (!!run.options['title']) ? run.options['title'].includes(filteringOptions.text): false;
      const descContainsText = (!filteringOptions.text.length) ? true : (!!run.options['description']) ? run.options['description'].includes(filteringOptions.text): false;

      const hasTags = (!filteringOptions.tags || filteringOptions.tags.length === 0) ? true :
        (!!run.options['tags']) ? filteringOptions.tags.every((searchedTag) => run.options['tags'].includes(searchedTag)): false;

      return isFavorite && (titleContainsText || descContainsText) && startedAfter && hasTags;
    });
  }

  private updateHistoryPane(historyRuns: DG.FuncCall[]) {
    this.historyList.updateList(historyRuns);
    this.historyFilter.updateTagList(historyRuns);
  };

  public updateRun(updatedRun: DG.FuncCall) {
    this.historyList.updateItem(updatedRun);
  }

  public addRun(newRun: DG.FuncCall) {
    this.historyList.addItem(newRun);

    this.allRuns.value.push(newRun);
  }

  public get history() {
    return this.allRuns.value;
  }

  constructor(
    private func: DG.Func,
  ) {
    super(ui.box(null, {style: {height: '100%'}}));

    this.root.appendChild(this.panel);

    const updateHistoryPaneSub = this.historyFilter.onFilteringChanged.subscribe(() => this.historyList.updateList(this.filteredHistoryRuns));

    const clickedSub = this.historyList.onClicked.subscribe((clickedCall) => this.onRunChosen.next(clickedCall.id));

    const comparisonSub = this.historyList.onComparisonCalled.subscribe((ids) => this.onComparison.next(ids));

    const allRunsSub = this.allRuns.subscribe(() => this.updateHistoryPane(this.filteredHistoryRuns));

    const allRunsFetch = this.allRunsFetch.subscribe(async () => {
      ui.setUpdateIndicator(this.root, true);
      const allRuns = (await historyUtils.pullRunsByName(this.func.name, [], {order: 'started'}, ['session.user', 'options'])).reverse();
      this.allRuns.next(allRuns);
      ui.setUpdateIndicator(this.root, false);
    });

    const onMetadataEdit = this.historyList.onMetadataEdit.subscribe((editedCall) => {
      {
        const run = this.allRuns.value.find((call) => call.id === editedCall.id);
        if (run) run.options = {...editedCall.options};
      }

      this.onRunEdited.next(editedCall);

      this.historyFilter.addTag(editedCall);
    });

    const onDeleteSub = this.historyList.onDelete.subscribe((deleteCall) => {
      const runIdx = this.allRuns.value.findIndex((call) => call.id === deleteCall.id);
      if (runIdx >=0) this.allRuns.value.splice(runIdx, 1);
      this.historyFilter.removeTag(deleteCall);

      this.afterRunDeleted.next(deleteCall);
    });

    this.subs.push(
      clickedSub,
      allRunsSub, allRunsFetch,
      onMetadataEdit,
      comparisonSub,
      onDeleteSub,
      updateHistoryPaneSub,
    );

    this.allRunsFetch.next();
  }
}
