/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import wu from 'wu';
import {BehaviorSubject, Subject, fromEvent, merge, of} from 'rxjs';
import {historyUtils} from '../../history-utils';
import '../css/history-panel.css';
import {getObservable, properUpdateIndicator} from '../../function-views/src/shared/utils';

class HistoricalRunsDelete extends DG.Dialog {
  public onFuncCallDelete = new Subject<Set<DG.FuncCall>>();

  constructor(funcCalls: Set<DG.FuncCall>) {
    const dlg = ui.dialog({title: `Delete ${funcCalls.size > 1 ? funcCalls.size: ''} ${funcCalls.size > 1 ? 'runs': 'run'}`});
    super(dlg.dart);

    this.add(ui.divText(`Deleted ${funcCalls.size > 1 ? funcCalls.size: ''} ${funcCalls.size > 1 ? 'runs': 'run'} ${funcCalls.size > 1 ? 'are': 'is'} impossible to restore. Are you sure?`));

    this.onOK(async () => this.onFuncCallDelete.next(funcCalls));
  }
}

export class HistoricalRunsList extends DG.Widget {
  public onComparisonCalled = new Subject<string[]>();

  public onMetadataEdit = new Subject<DG.FuncCall>();
  public onDelete = new Subject<DG.FuncCall>();
  public onClicked = new Subject<DG.FuncCall>();

  private cards = this.initialFuncCalls.map((funcCall) => new HistoricalRunCard(funcCall));

  private onSelectedRunsChanged = new BehaviorSubject<Set<DG.FuncCall>>(new Set());
  private onRunsChanged = new BehaviorSubject<DG.FuncCall[]>(this.initialFuncCalls);

  private actions = ui.box() as HTMLElement;

  public get selected() {
    return this.onSelectedRunsChanged.value;
  }

  constructor(private initialFuncCalls: DG.FuncCall[], private options?: {
    fallbackText?: string,
  }) {
    super(ui.div());

    const listChangedSub = this.onRunsChanged.subscribe((runs) => {
      this.refresh(runs);
      this.onSelectedRunsChanged.next(new Set([]));
    });
    const selectionChangedSub = this.onSelectedRunsChanged.subscribe((selectedRuns) => {
      this.cards.forEach((card) => {
        if (selectedRuns.has(card.funcCall))
          card.setSelected();
        else
          card.setUnselected();
      });
      this.refreshActions();
    });
    this.subs.push(listChangedSub, selectionChangedSub);
  }

  deleteItem(id: string) {
    const corrCards = this.cards.splice(this.cards.findIndex((card) => card.funcCall.id === id)!, 1);
    // We delete single item
    corrCards[0].deleteCard(false);
  }

  updateItem(updatedFuncCall: DG.FuncCall) {
    const corrCard = this.cards.find((card) => card.funcCall.id === updatedFuncCall.id);
    corrCard?.updateCard(updatedFuncCall);
  }

  addItem(newFuncCall: DG.FuncCall) {
    this.onRunsChanged.next([newFuncCall, ...this.onRunsChanged.value]);
  }

  updateList(newFuncCalls: DG.FuncCall[]) {
    this.onRunsChanged.next(newFuncCalls);
  }

  private get getCompareIcon() {
    const currentSelectedSet = this.onSelectedRunsChanged.value;

    const t = ui.iconFA('exchange', async () => {
      this.onComparisonCalled.next([...wu(currentSelectedSet.keys()).map((selected) => selected.id)]);
    }, 'Compare selected runs');
    t.style.margin = '5px';
    if (currentSelectedSet.size < 2)
      t.classList.add('hp-disabled');
    return t;
  }

  private get getTrashIcon() {
    const t = ui.iconFA('trash-alt', () => {
      const deleteDialog = new HistoricalRunsDelete(this.onSelectedRunsChanged.value);

      const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async (setToDelete) => {
        await Promise.all(
          wu(setToDelete.values()).map(async (funcCall) => {
            this.deleteItem(funcCall.id);
            await historyUtils.deleteRun(funcCall);

            return Promise.resolve();
          }));

        onDeleteSub.unsubscribe();
      });

      deleteDialog.show({center: true, width: 500});
    }, 'Delete selected runs');

    t.style.margin = '5px';
    if (this.onSelectedRunsChanged.value.size === 0)
      t.classList.add('hp-disabled');
    return t;
  }

  private get selectAllIcon() {
    const t = ui.iconFA('square', () => this.onSelectedRunsChanged.next(new Set(this.onRunsChanged.value)), 'Select all');
    t.style.margin = '5px';
    return t;
  }

  private get unselectAllIcon() {
    const fullListCount = this.onRunsChanged.value.length;
    const currentSelectedSet = this.onSelectedRunsChanged.value;
    const iconType = currentSelectedSet.size === fullListCount? 'check-square': 'minus-square';
    const t = ui.iconFA(iconType, () => this.onSelectedRunsChanged.next(new Set()), 'Unselect all');
    t.style.margin = '5px';
    return t;
  }

  private refresh(newFuncCalls: DG.FuncCall[]) {
    ui.empty(this.root);
    this.cards = newFuncCalls.map((funcCall) => new HistoricalRunCard(funcCall));

    if (newFuncCalls.length > 0) {
      this.root.appendChild(ui.divV([
        this.actions,
        ui.element('div', 'splitbar-horizontal'),
        ...this.cards.map((card) => card.root),
      ]));
    } else {
      this.root.appendChild(ui.divV([
        ui.element('div', 'splitbar-horizontal'),
        ui.divText(this.options?.fallbackText ?? 'No historical runs found', 'hp-no-elements-label'),
      ]));
    }

    const anyClickSub = merge(...this.cards.map((card) => card.onCardClicked))
      .subscribe((clickedCall) => this.onClicked.next(clickedCall));

    const anySelectSub = merge(...this.cards.map((card) => card.onSelect))
      .subscribe((selectedCall) => this.onSelectedRunsChanged.next(new Set([
        selectedCall, ...this.onSelectedRunsChanged.value,
      ])));

    const anyUnselectSub = merge(...this.cards.map((card) => card.onUnselect))
      .subscribe((selectedCall) => {
        this.onSelectedRunsChanged.value.delete(selectedCall);
        this.onSelectedRunsChanged.next(this.onSelectedRunsChanged.value);
      });

    const anyDeleteSub = merge(...this.cards.map((card) => card.onDelete))
      .subscribe((callToDelete) => {
        this.onSelectedRunsChanged.value.delete(callToDelete);
        this.onSelectedRunsChanged.next(this.onSelectedRunsChanged.value);

        this.onDelete.next(callToDelete);
      });

    const anyMetadataEditSub = merge(...this.cards.map((card) => card.onMetadataEdit))
      .subscribe((editedCall) => {
        this.onMetadataEdit.next(editedCall);
      });
    this.subs.push(anyClickSub, anySelectSub, anyUnselectSub, anyDeleteSub, anyMetadataEditSub);
  }

  private refreshActions() {
    const newActions = this.getActions;
    this.actions.replaceWith(newActions);
    this.actions = newActions;
  }

  private get getActions() {
    const currentSelectedSet = this.onSelectedRunsChanged.value;
    const actionsSection = ui.divH([
      ui.span([`Selected: ${currentSelectedSet.size}`], {style: {'align-self': 'center'}}),
      ui.divH([
        this.getCompareIcon,
        this.getTrashIcon,
        currentSelectedSet.size === 0 ? this.selectAllIcon: this.unselectAllIcon,
      ]),
    ], {style: {
      'justify-content': 'space-between',
      'padding': '0 12px',
    }}) as HTMLElement;

    return actionsSection;
  }
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
    dummyInput.root.style.display = 'none';

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

      if (tags!.length === 0)
        dummyInput.root.style.display = 'none';
      else
        dummyInput.root.style.removeProperty('display');
    });

    const updateTagsSub = this.onFuncCallListChanged.subscribe((newFuncCalls) => {
      properUpdateIndicator(tagInput.root, true);
      const newTags = newFuncCalls
        .map((run) => run.options['tags'] as string[])
        .reduce((acc, runTags) => {
          if (!!runTags) runTags.forEach((runTag) => {if (!acc.includes(runTag)) acc.push(runTag);});
          return acc;
        }, [] as string[]);

      tagInput.items = ['Choose tag to filter', ...newTags];
      properUpdateIndicator(tagInput.root, false);
    });

    this.subs.push(filterChangeSub, addTagSub, updateTagsSub);

    const form = ui.divV([
      favoritesInput,
      textInput,
      dateInput,
      tagInput,
      dummyInput,
    ], 'ui-form ui-form-wide ui-form-left');

    this.root.appendChild(form);
  }

  addTag(funccall: DG.FuncCall) {
    const newFunccalls = this.onFuncCallListChanged.value.filter((fc) => fc.id !== funccall.id);
    this.onFuncCallListChanged.next([...newFunccalls, funccall]);
  }

  removeTag(funccall: DG.FuncCall) {
    const newFunccalls = this.onFuncCallListChanged.value.filter((fc) => fc.id !== funccall.id);
    this.onFuncCallListChanged.next(newFunccalls);
  }

  updateTagList(funccalls: DG.FuncCall[]) {
    this.onFuncCallListChanged.next(funccalls);
  }
}

type EditOptions = {
  title: string | null,
  description: string | null,
  tags: string[],
  favorite: 'favorited' | 'unfavorited' | 'same'
}

export class HistoricalRunEdit extends DG.Dialog {
  public onMetadataEdit = new Subject<EditOptions>();

  constructor(funcCall: DG.FuncCall) {
    const dlg = ui.dialog({title: 'Edit run metadata'});

    super(dlg.dart);

    let title = funcCall.options['title'] ?? '';
    let description = funcCall.options['description'] ?? '';
    let isFavorite = funcCall.options['isFavorite'] ?? false;
    const titleInput = ui.stringInput('Title', title, (s: string) => title = s);

    const dummyInput = ui.stringInput(' ', '');
    const tagsLine = DG.TagEditor.create();
    (funcCall.options['tags'] ?? []).forEach((tag: string) => {
      tagsLine.addTag(tag);
    });
    dummyInput.input.replaceWith(tagsLine.root);

    const addNewTag = () => {
      if (tagInput.value === '' ||
        // @ts-ignore
        tagsLine.tags.includes(tagInput.value)
      )
        return;

      tagsLine.addTag(tagInput.value);
      tagInput.value = '';
    };
    const tagInput = ui.stringInput('Tag', '').addOptions(ui.iconFA('plus', addNewTag, 'Add this tag'));

    const enterSub = fromEvent<KeyboardEvent>(tagInput.input, 'onkeydown').subscribe((ev) => {
      if (ev.key == 'Enter') {
        ev.stopPropagation();
        addNewTag();
      }
    });
    this.sub(enterSub);

    this.add(ui.form([
      titleInput,
      ui.stringInput('Description', description, (s: string) => description = s),
      ui.boolInput('Favorites', isFavorite, (b: boolean) => isFavorite = b),
      tagInput,
      dummyInput,
    ]));

    this.addButton('Save', async () => {
      if (tagInput.value !== '') {
        grok.shell.info(`Dialog has unsaved tags: ${tagInput.value}`);
        return;
      }

      let favorite = 'same';
      if (isFavorite && !funcCall.options['isFavorite']) favorite = 'favorited';
      if (!isFavorite && funcCall.options['isFavorite']) favorite = 'unfavorited';

      const editOptions = {
        title: (title !== '') ? title : null,
        description: (description !== '') ? description : null,
        tags: tagsLine.tags.map((el) => el as any),
        favorite,
      } as EditOptions;

      this.onMetadataEdit.next(editOptions);
      this.close();
    });
  }
}

class HistoricalRunCard extends DG.Widget {
  public onCardClicked = new Subject<DG.FuncCall>();
  public onSelect = new Subject<DG.FuncCall>();
  public onUnselect = new Subject<DG.FuncCall>();
  public onMetadataEdit = new Subject<DG.FuncCall>();
  public onDelete = new Subject<DG.FuncCall>();

  private onFuncCallChanged = new BehaviorSubject<DG.FuncCall>(this.initialFuncCall);

  public get funcCall() {
    return this.onFuncCallChanged.value;
  }

  private addToSelected = ui.iconFA('square', (ev) => {
    ev.stopPropagation();
    this.onSelect.next(this.funcCall);
    this.setSelected();
  }, 'Select this run');

  private removeFromSelected = ui.iconFA('check-square', (ev) => {
    ev.stopPropagation();
    this.onUnselect.next(this.funcCall);
    this.setUnselected();
  }, 'Unselect this run');

  private addToFavorites = ui.iconFA('star', async (ev) => {
    ev.stopPropagation();
    this.setFavorited();
  }, 'Add to favorites');

  private unfavoriteIcon = ui.iconFA('star', async (ev) => {
    ev.stopPropagation();
    this.setUnfavorited();
  }, 'Unfavorite the run');

  constructor(
    private initialFuncCall: DG.FuncCall,
    private options = {
      showEdit: true,
      showDelete: true,
      showFavorite: true,
      showSelect: true,
      showAuthorIcon: true,
    },
  ) {
    super(ui.div());

    const onChangedSub = this.onFuncCallChanged.subscribe((funccall) => this.refresh(funccall));

    const clickSub = fromEvent(this.root, 'click').subscribe(async () => {
      this.onCardClicked.next(initialFuncCall);

      this.root.classList.add('clicked');
    });
    this.subs.push(clickSub, onChangedSub);
  }

  async setFavorited() {
    properUpdateIndicator(this.root, true);
    return historyUtils.loadRun(this.funcCall.id, false).then((fullCall) => {
      fullCall.options['isFavorite'] = true;
      return historyUtils.saveRun(fullCall);
    }).then((fullCall) => {
      this.addToFavorites.style.display = 'none';
      this.unfavoriteIcon.style.removeProperty('display');

      this.onMetadataEdit.next(fullCall);
    }).catch((err) => {
      grok.shell.error(err);
    }).finally(() => {
      properUpdateIndicator(this.root, false);
    });
  }

  async setUnfavorited() {
    properUpdateIndicator(this.root, true);
    return historyUtils.loadRun(this.funcCall.id).then((fullCall) => {
      fullCall.options['isFavorite'] = false;
      return historyUtils.saveRun(fullCall);
    }).then((fullCall) => {
      this.unfavoriteIcon.style.display = 'none';
      this.addToFavorites.style.removeProperty('display');

      this.onMetadataEdit.next(fullCall);
    }).catch((err) => {
      grok.shell.error(err);
    }).finally(() => {
      properUpdateIndicator(this.root, false);
    });
  }

  setSelected() {
    this.addToSelected.style.display = 'none';
    this.removeFromSelected.style.removeProperty('display');
  }

  setUnselected() {
    this.removeFromSelected.style.display = 'none';
    this.addToSelected.style.removeProperty('display');
  }

  updateCard(updatedFunccall: DG.FuncCall) {
    this.onFuncCallChanged.next(updatedFunccall);
    this.onMetadataEdit.next(updatedFunccall);
  }

  async deleteCard(showDialog = true) {
    if (showDialog) {
      const deleteDialog = new HistoricalRunsDelete(new Set([this.funcCall]));

      const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async (setToDelete) => {
        properUpdateIndicator(this.root, true);
        Promise.all(
          wu(setToDelete.values()).map(async () => {
            return historyUtils.deleteRun(this.funcCall).then(() => {
              this.onDelete.next(this.funcCall);
            });
          }))
          .then(() => {
            ui.empty(this.root);
            onDeleteSub.unsubscribe();
          })
          .catch((err) => {
            grok.shell.error(err);
          }).finally(() => {
            properUpdateIndicator(this.root, false);
          });
      });
      deleteDialog.show({center: true, width: 500});
    } else {
      historyUtils.deleteRun(this.funcCall).then(() => {
        ui.empty(this.root);

        this.onDelete.next(this.funcCall);
      });
    }

    return this.onDelete.toPromise();
  }

  async editCard() {
    const editDialog = new HistoricalRunEdit(this.funcCall);

    const onEditSub = editDialog.onMetadataEdit.subscribe(async (editOptions) => {
      properUpdateIndicator(this.root, true);
      return historyUtils.loadRun(this.funcCall.id, false)
        .then((fullCall) => {
          if (editOptions.title) fullCall.options['title'] = editOptions.title;
          if (editOptions.description) fullCall.options['description'] = editOptions.description;
          if (editOptions.tags) fullCall.options['tags'] = editOptions.tags;
          if (editOptions.favorite !== 'same') fullCall.options['isFavorite'] = (editOptions.favorite === 'favorited');

          return historyUtils.saveRun(fullCall);
        })
        .then((fullCall) => {
          this.onFuncCallChanged.next(fullCall);
          this.onMetadataEdit.next(fullCall);

          onEditSub.unsubscribe();
        })
        .catch((err) => {
          grok.shell.error(err);
        }).finally(() => {
          properUpdateIndicator(this.root, false);
        });
    });
    editDialog.show({center: true, width: 500});

    return this.onMetadataEdit.toPromise();
  }

  private refresh(funcCall: DG.FuncCall) {
    const options = this.options;

    const icon = funcCall.author.picture as HTMLElement;
    icon.style.width = '25px';
    icon.style.height = '25px';
    icon.style.fontSize = '20px';
    const cardLabel = ui.label(funcCall.options['title'] ?? funcCall.author.friendlyName, {style: {'color': 'var(--blue-1)'}});
    ui.bind(funcCall.author, icon);

    const editIcon = ui.iconFA('edit', (ev) => {
      ev.stopPropagation();
      this.editCard().then(() =>this.onMetadataEdit.next(funcCall));
    }, 'Edit run metadata');
    editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
      ev.stopPropagation();
      this.deleteCard().then(() => this.onDelete.next(funcCall));
    }, 'Delete run');
    deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    this.addToFavorites.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.unfavoriteIcon.classList.add('fas', 'hp-funccall-card-icon');

    if (funcCall.options['isFavorite']) {
      this.addToFavorites.style.display = 'none';
      this.unfavoriteIcon.style.removeProperty('display');
    } else {
      this.unfavoriteIcon.style.display = 'none';
      this.addToFavorites.style.removeProperty('display');
    }

    this.addToSelected.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.removeFromSelected.classList.add('hp-funccall-card-icon');
    this.removeFromSelected.style.display = 'none';

    const dateStarted = new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});

    const card = ui.divH([
      ui.divH([
        ...(options.showAuthorIcon) ? [icon]: [],
        ui.divV([
          cardLabel,
          ui.span([dateStarted]),
          ...(funcCall.options['description']) ? [ui.divText(funcCall.options['description'], 'description')]: [],
          ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
            [ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')))]:[],
        ], 'hp-card-content'),
      ]),
      ui.divH([
        ...(options.showFavorite) ? [this.unfavoriteIcon, this.addToFavorites]: [],
        ...(options.showEdit) ? [editIcon]: [],
        ...(options.showDelete) ? [deleteIcon]: [],
        ...(options.showSelect) ? [this.addToSelected, this.removeFromSelected]: [],
      ]),
    ], 'hp-funccall-card');

    ui.tooltip.bind(card, () => ui.tableFromMap({
      Author: grok.shell.user.toMarkup(),
      Date: dateStarted,
      ...(funcCall.options['title']) ? {'Title': funcCall.options['title']}:{},
      ...(funcCall.options['description']) ? {'Description': funcCall.options['description']}:{},
      ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
        {'Tags': ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor')}:{},
    }));

    ui.empty(this.root);
    this.root.appendChild(card);
  }
}

type FilterOptions = {
  isOnlyFavorites: boolean,
  text: string,
  author?: DG.User,
  startedAfter?: dayjs.Dayjs,
  tags?: string[]
}

export class HistoryPanel extends DG.Widget {
  public addRun(newRun: DG.FuncCall) {
    this.historyList.addItem(newRun);
  }

  // Emitted when FuncCall should is chosen. Contains FuncCall ID
  public onRunChosen = new Subject<string>();

  // Emitted when FuncCalls are called for comparison. Contains FuncCalls' IDs
  public onComparison = new Subject<string[]>();

  // Emitted when FuncCall is edited
  public onRunEdited = new Subject<DG.FuncCall>();

  // Emitted when FuncCall is deleted
  public afterRunDeleted = new Subject<DG.FuncCall>();

  public allRunsFetch = new Subject<true>();

  historyFilter = new HistoryFilter();

  historyList = new HistoricalRunsList([], {fallbackText: 'No runs are found in history'});

  panel = ui.divV([
    this.historyFilter.root,
    this.historyList.root,
  ]);

  updateHistoryPane(myRuns: DG.FuncCall[]) {
    this.historyList.updateList(myRuns);
    this.historyFilter.updateTagList(myRuns);
  };

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
      this.allRuns.value.splice(runIdx, 1);
      this.historyFilter.removeTag(deleteCall);
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
