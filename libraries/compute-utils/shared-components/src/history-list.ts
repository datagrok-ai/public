import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import {Subject, BehaviorSubject, merge, fromEvent} from 'rxjs';
import {historyUtils} from '../../history-utils';
import {properUpdateIndicator} from '../../function-views/src/shared/utils';
import {HistoricalRunsDelete, HistoricalRunEdit} from './history-dialogs';
import {EXPERIMENTAL_TAG, storageName} from '../../shared-utils/consts';
import {FilterOptions} from './history-panel';

class HistoricalRunCard extends DG.Widget {
  private _onFuncCallChanged = new BehaviorSubject<DG.FuncCall>(this.initialFunccall);
  private _onFavoriteChanged = new BehaviorSubject<boolean>(this.initialIsFavorite);
  private favStorageName = `${storageName}_${this.funcCall.func.name}_Fav`;

  private _onClicked = new Subject<DG.FuncCall>();
  private _onSelect = new Subject<DG.FuncCall>();
  private _onUnselect = new Subject<DG.FuncCall>();
  private _onMetadataEdit = new Subject<DG.FuncCall>();
  private _onDelete = new Subject<DG.FuncCall>();

  private addToSelected = ui.iconFA('square', (ev) => {
    ev.stopPropagation();
    this.redrawIsSelected(true);
    this._onSelect.next(this.funcCall);
  }, 'Select this run');

  private removeFromSelected = ui.iconFA('check-square', (ev) => {
    ev.stopPropagation();
    this.redrawIsSelected(false);
    this._onUnselect.next(this.funcCall);
  }, 'Unselect this run');

  private addToFavorites = ui.iconFA('star', async (ev) => {
    ev.stopPropagation();
    this.saveIsFavorite(true);
  }, 'Add to favorites');

  private unfavoriteIcon = ui.iconFA('star', async (ev) => {
    ev.stopPropagation();
    this.saveIsFavorite(false);
  }, 'Unfavorite the run');

  public onClicked = this._onClicked.asObservable();
  public onSelect = this._onSelect.asObservable();
  public onUnselect = this._onUnselect.asObservable();
  public onFavoriteChanged = this._onFavoriteChanged.asObservable();
  public onMetadataEdit = this._onMetadataEdit.asObservable();
  public onDelete = this._onDelete.asObservable();

  public get funcCall() {
    return this._onFuncCallChanged.value;
  }

  public get onFuncCallChanged() {
    return this._onFuncCallChanged.asObservable();
  }

  public get isFavorite() {
    return this._onFavoriteChanged.value;
  }

  constructor(
    private readonly initialFunccall: DG.FuncCall,
    private readonly initialIsFavorite: boolean,
    private readonly options?: {
      showEdit?: boolean,
      showDelete?: boolean,
      showFavorite?: boolean,
      showAuthorIcon?: boolean,
    },
  ) {
    super(ui.div(ui.divText('No run is loaded', 'hp-no-elements-label')));

    const onChangedSub = merge(this._onFuncCallChanged, this.onMetadataEdit)
      .subscribe((funccall) => this.redraw(funccall, this.isFavorite));

    const clickSub = fromEvent(this.root, 'click').subscribe(() => this._onClicked.next(this.funcCall));
    this.subs.push(clickSub, onChangedSub);
  }

  private async saveIsFavorite(isFavorite: boolean) {
    properUpdateIndicator(this.root, true);
    if (isFavorite) {
      return grok.dapi.userDataStorage.postValue(this.favStorageName, this.funcCall.id, '')
        .then(() => {
          this.redrawIsFavorite(true);
          this._onFavoriteChanged.next(true);
        })
        .finally(() => properUpdateIndicator(this.root, false));
    } else {
      return grok.dapi.userDataStorage.remove(this.favStorageName, this.funcCall.id)
        .then(() => {
          this.redrawIsFavorite(false);
          this._onFavoriteChanged.next(false);
        })
        .finally(() => properUpdateIndicator(this.root, false));
    }
  }

  redrawIsFavorite(isFavorite: boolean) {
    if (isFavorite) {
      ui.setDisplay(this.addToFavorites, false);
      ui.setDisplay(this.unfavoriteIcon, true);
    } else {
      ui.setDisplay(this.unfavoriteIcon, false);
      ui.setDisplay(this.addToFavorites, true);
    }
  }

  redrawIsSelected(isSelected: boolean) {
    if (isSelected) {
      ui.setDisplay(this.addToSelected, false);
      ui.setDisplay(this.removeFromSelected, true);
    } else {
      ui.setDisplay(this.removeFromSelected, false);
      ui.setDisplay(this.addToSelected, true);
    }
  }

  redrawRun(updatedRun: DG.FuncCall) {
    this._onFuncCallChanged.next(updatedRun);
    this._onMetadataEdit.next(updatedRun);
  }

  async showDeleteDialog() {
    const deleteDialog = new HistoricalRunsDelete(new Set([this.funcCall]));

    const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(() => {
      this._onDelete.next(this.funcCall);
      onDeleteSub.unsubscribe();
    });
    deleteDialog.show({center: true, width: 500});
  }

  async showEditDialog() {
    const editDialog = new HistoricalRunEdit(this.funcCall, this.isFavorite);

    const onEditSub = editDialog.onMetadataEdit.subscribe(async (editOptions) => {
      properUpdateIndicator(this.root, true);

      return ((editOptions.favorite !== 'same') ?
        this.saveIsFavorite((editOptions.favorite === 'favorited')) :
        Promise.resolve())
        .then(() => historyUtils.loadRun(this.funcCall.id, false))
        .then((fullCall) => {
          if (editOptions.title) fullCall.options['title'] = editOptions.title;
          if (editOptions.description) fullCall.options['description'] = editOptions.description;
          if (editOptions.tags) fullCall.options['tags'] = editOptions.tags;

          return Promise.all([historyUtils.saveRun(fullCall), fullCall]);
        })
        .then(([, fullCall]) => {
          this._onMetadataEdit.next(fullCall);
          this._onFuncCallChanged.next(fullCall);

          onEditSub.unsubscribe();
        })
        .catch((err) => {
          grok.shell.error(err);
        }).finally(() => {
          properUpdateIndicator(this.root, false);
        });
    });
    editDialog.show({center: true, width: 500});
  }

  private redraw(run: DG.FuncCall, isFavorite: boolean) {
    const options = this.options;

    let authorIcon = null;
    if (this.options?.showAuthorIcon) {
      authorIcon = run.author.picture as HTMLElement;
      authorIcon.style.width = '25px';
      authorIcon.style.height = '25px';
      authorIcon.style.fontSize = '20px';

      ui.bind(run.author, authorIcon);
    }

    const experimentalTag = ui.iconFA('flask', null, 'Experimental run');
    experimentalTag.classList.add('fad', 'fa-sm');
    experimentalTag.classList.remove('fal');
    experimentalTag.style.marginLeft = '3px';
    const immutableTags = run.options['immutable_tags'] as string[] | undefined;
    const cardLabel = ui.span([
      ui.label(
        run.options['title'] ??
        run.author?.friendlyName ??
        grok.shell.user.friendlyName, {style: {'color': 'var(--blue-1)'}},
      ),
      ...(immutableTags && immutableTags.includes(EXPERIMENTAL_TAG)) ?
        [experimentalTag]:[],
    ]);

    const editIcon = ui.iconFA('edit', (ev) => {
      ev.stopPropagation();
      this.showEditDialog();
    }, 'Edit run metadata');
    editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
      ev.stopPropagation();
      this.showDeleteDialog();
    }, 'Delete run');
    deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    this.addToFavorites.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.unfavoriteIcon.classList.add('fas', 'hp-funccall-card-icon');

    if (isFavorite) {
      ui.setDisplay(this.addToFavorites, false);
      ui.setDisplay(this.unfavoriteIcon, true);
    } else {
      ui.setDisplay(this.unfavoriteIcon, false);
      ui.setDisplay(this.addToFavorites, true);
    }

    this.addToSelected.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.removeFromSelected.classList.add('hp-funccall-card-icon');
    ui.setDisplay(this.removeFromSelected, false);

    const dateStarted = new Date(run.started.toString())
      .toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});

    const card = ui.divH([
      ui.divH([
        ...(options?.showAuthorIcon) ? [authorIcon!]: [],
        ui.divV([
          cardLabel,
          ui.span([dateStarted]),
          ...(run.options['description']) ? [ui.divText(run.options['description'], 'description')]: [],
          ...(run.options['tags'] && run.options['tags'].length > 0) ?
            [ui.div(run.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')))]:[],
        ], 'hp-card-content'),
      ]),
      ui.divH([
        ...(options?.showFavorite) ? [this.unfavoriteIcon, this.addToFavorites]: [],
        ...(options?.showEdit) ? [editIcon]: [],
        ...(options?.showDelete) ? [deleteIcon]: [],
        this.addToSelected, this.removeFromSelected,
      ]),
    ], 'hp-funccall-card');

    ui.tooltip.bind(card, () => ui.tableFromMap({
      Author: grok.shell.user.toMarkup(),
      Date: dateStarted,
      ...(run.options['title']) ? {'Title': run.options['title']}:{},
      ...(run.options['description']) ? {'Description': run.options['description']}:{},
      ...(run.options['tags'] && run.options['tags'].length > 0) ?
        {'Tags': ui.div(run.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor')}:{},
      ...(run.options['immutable_tags'] && run.options['immutable_tags'].length > 0) ?
        {'Extra': ui.div(run.options['immutable_tags']
          .map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor',
        )}:{},
    }));

    ui.empty(this.root);
    this.root.appendChild(card);
  }
}

export class HistoricalRunsList extends DG.Widget {
  private storageName(runs: DG.FuncCall[]) {
    return `${storageName}_${runs[0].func.name}_Fav`;
  }

  public onComparisonCalled = new Subject<string[]>();
  public onMetadataEdit = new Subject<DG.FuncCall>();
  public onDelete = new Subject<DG.FuncCall>();
  public onClicked = new Subject<DG.FuncCall>();

  private _onSelectedChanged = new BehaviorSubject<Set<DG.FuncCall>>(new Set([]));
  private onFilteringChanged = new BehaviorSubject<FilterOptions>({isOnlyFavorites: false, text: ''});
  private onRunsChanged = new BehaviorSubject<DG.FuncCall[]>(this.initialRuns);
  private onFilteredRunsChanged = new BehaviorSubject<DG.FuncCall[]>(this.onRunsChanged.value);

  private cards: HistoricalRunCard[] = [];
  private total = ui.span([`Selected: ${this._onSelectedChanged.value.size}`], {style: {'align-self': 'center'}});

  public get runs() {
    return this.onRunsChanged.value;
  }

  private isOKforFilter(card: HistoricalRunCard) {
    const filteringOptions = this.onFilteringChanged.value;
    const run = card.funcCall;

    const isFavorite = filteringOptions.isOnlyFavorites ? card.isFavorite : true;
    const startedAfter = (filteringOptions.startedAfter ? run.started > filteringOptions.startedAfter : true);
    const titleContainsText = (!filteringOptions.text.length) ? true :
      (!!run.options['title']) ? run.options['title'].includes(filteringOptions.text): false;
    const descContainsText = (!filteringOptions.text.length) ? true :
      (!!run.options['description']) ? run.options['description'].includes(filteringOptions.text): false;

    const hasTags = (!filteringOptions.tags || filteringOptions.tags.length === 0) ? true :
      (!!run.options['tags']) ?
        filteringOptions.tags.every((searchedTag) => run.options['tags'].includes(searchedTag)): false;

    return isFavorite && (titleContainsText || descContainsText) && startedAfter && hasTags;
  }

  public filter() {
    this.cards.forEach((card) => ui.setDisplay(card.root, this.onFilteredRunsChanged.value.includes(card.funcCall)));
  }

  public get selected() {
    return this._onSelectedChanged.value;
  }

  public get onSelectedChanged() {
    return this._onSelectedChanged.asObservable();
  }

  constructor(
    private readonly initialRuns: DG.FuncCall[],
    private readonly options: {
      fallbackText?: string,
      showEdit?: boolean,
      showDelete?: boolean,
      showFavorite?: boolean,
      showAuthorIcon?: boolean,
      showCompare?: boolean,
    } = {
      showEdit: true,
      showDelete: true,
      showFavorite: true,
      showAuthorIcon: true,
      showCompare: true,
    }) {
    super(ui.div([
      ui.divV([
        ui.element('div', 'splitbar-horizontal'),
        ui.divText('No historical runs were loaded', 'hp-no-elements-label'),
      ]),
    ], {style: {'height': '100%', 'overflow-y': 'hidden'}}));

    const selectionChangedSub = this._onSelectedChanged.subscribe((selectedRuns) => {
      this.cards.forEach((card) => {
        card.redrawIsSelected(selectedRuns.has(card.funcCall));
      });
      this.redrawSelectionState();
    });
    const listChangedSub = this.onRunsChanged.subscribe(async (runs) => {
      if (runs.length === 0)
        this.redraw(runs, []);
      else {
        grok.dapi.userDataStorage
          .get(this.storageName(runs))
          .then((data: Record<string, string>) => {
            this.redraw(runs, Object.keys(data));
          })
          .catch(() => {
            grok.shell.error(`Favorites storage is not available. Try again later.`);
            this.redraw(runs, []);
          })
          .finally(() => {
            this._onSelectedChanged.next(new Set([]));
          });
      }
    });

    const filteringChangedSub = this.onFilteringChanged.subscribe(() => {
      this.onFilteredRunsChanged.next(
        this.cards.filter((card) => this.isOKforFilter(card)).map((card) => card.funcCall),
      );
      this._onSelectedChanged.next(new Set([]));
      this.filter();
    });

    this.subs.push(listChangedSub, selectionChangedSub, filteringChangedSub);
  }

  async deleteRun(id: string) {
    const corrCards = this.cards.splice(this.cards.findIndex((card) => card.funcCall.id === id), 1);
    const deletedCard = corrCards[0];
    properUpdateIndicator(deletedCard.root, true);
    return historyUtils.deleteRun(deletedCard.funcCall).then(() => {
      ui.empty(deletedCard.root);
      const idx = this.runs.findIndex((run) => run.id === id);
      if (idx >= 0)
        this.runs.splice(idx, 1);
    }).catch((e) => {
      grok.shell.error(e);
    }).finally(() => {
      properUpdateIndicator(deletedCard.root, false);
    });
  }

  editRun(editedCall: DG.FuncCall) {
    const corrCard = this.cards.find((card) => card.funcCall.id === editedCall.id);
    corrCard?.showEditDialog();
  }

  updateRun(updatedRun: DG.FuncCall) {
    const corrCard = this.cards.find((card) => card.funcCall.id === updatedRun.id);
    corrCard?.redrawRun(updatedRun);
  }

  addRun(newRun: DG.FuncCall) {
    this.onRunsChanged.next([newRun, ...this.onRunsChanged.value]);
  }

  updateRuns(newRuns: DG.FuncCall[]) {
    this.onRunsChanged.next(newRuns);
  }

  setFiltering(filteringOptions: FilterOptions) {
    this.onFilteringChanged.next(filteringOptions);
  }

  resetFiltering() {
    this.onFilteringChanged.next({isOnlyFavorites: false, text: ''});
  }

  private getTrashIcon() {
    const t = ui.iconFA('trash-alt', () => {
      const deleteDialog = new HistoricalRunsDelete(this._onSelectedChanged.value);

      const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async (setToDelete) => {
        await Promise.all(
          wu(setToDelete.values()).map(async (funcCall) => {
            await this.deleteRun(funcCall.id);

            this._onSelectedChanged.value.delete(funcCall);
            this._onSelectedChanged.next(this._onSelectedChanged.value);

            return Promise.resolve();
          }));

        onDeleteSub.unsubscribe();
      });

      deleteDialog.show({center: true, width: 500});
    }, 'Delete selected runs');

    t.style.margin = '5px';
    if (this._onSelectedChanged.value.size === 0)
      t.classList.add('hp-disabled');
    return t;
  }
  private trashIcon = this.getTrashIcon();

  private getSelectAllIcon() {
    const t = ui.iconFA('square',
      () => this._onSelectedChanged.next(new Set(this.onFilteredRunsChanged.value)),
      'Select all',
    );
    t.style.margin = '5px';
    return t;
  }
  private selectAllIcon = this.getSelectAllIcon();

  private getUnselectOnPartialIcon() {
    const t = ui.iconFA('minus-square', () => this._onSelectedChanged.next(new Set()), 'Unselect all');
    t.style.margin = '5px';
    return t;
  }
  private unselectOnPartialIcon = this.getUnselectOnPartialIcon();

  private getUnselectAllIcon() {
    const t = ui.iconFA('check-square', () => this._onSelectedChanged.next(new Set()), 'Unselect all');
    t.style.margin = '5px';
    return t;
  }
  private unselectAllIcon = this.getUnselectAllIcon();

  private getCompareIcon() {
    const t = ui.iconFA('exchange', async () => {
      this.onComparisonCalled.next([...wu(this._onSelectedChanged.value.keys()).map((selected) => selected.id)]);
    }, 'Compare selected runs');
    t.style.margin = '5px';
    return t;
  }
  private compareIcon = this.getCompareIcon();

  private redraw(newFuncCalls: DG.FuncCall[], favorites: string[]) {
    ui.empty(this.root);
    this.cards = newFuncCalls.map((funcCall) => new HistoricalRunCard(
      funcCall,
      favorites.includes(funcCall.id),
      this.options),
    );

    if (newFuncCalls.length > 0) {
      this.root.appendChild(ui.divV([
        ui.divH([
          this.total,
          ui.divH([
            ...this.options.showCompare ? [this.compareIcon]: [],
            this.trashIcon,
            this.selectAllIcon, this.unselectOnPartialIcon, this.unselectAllIcon,
          ]),
        ], {style: {
          'justify-content': 'space-between',
          'padding': '0px 10px 0px 9px',
        }}),
        ui.element('div', 'splitbar-horizontal'),
        ui.divV(this.cards.map((card) => card.root), {style: {'overflow-y': 'scroll'}}),
      ], {style: {'overflow-y': 'hidden', 'height': '100%'}}));
    } else {
      this.root.appendChild(ui.divV([
        ui.element('div', 'splitbar-horizontal'),
        ui.divText(this.options?.fallbackText ?? 'No historical runs found', 'hp-no-elements-label'),
      ]));
    }

    const anyClickSub = merge(...this.cards.map((card) => card.onClicked))
      .subscribe((clickedCall) => this.onClicked.next(clickedCall));

    const anySelectSub = merge(...this.cards.map((card) => card.onSelect))
      .subscribe((selectedCall) => this._onSelectedChanged.next(new Set([
        selectedCall, ...this._onSelectedChanged.value,
      ])));

    const anyUnselectSub = merge(...this.cards.map((card) => card.onUnselect))
      .subscribe((selectedCall) => {
        this._onSelectedChanged.value.delete(selectedCall);
        this._onSelectedChanged.next(this._onSelectedChanged.value);
      });

    const anyDeleteSub = merge(...this.cards.map((card) => card.onDelete))
      .subscribe((callToDelete) => {
        this.deleteRun(callToDelete.id).then(() => {
          this._onSelectedChanged.value.delete(callToDelete);
          this._onSelectedChanged.next(this._onSelectedChanged.value);

          const runIdx = this.onRunsChanged.value.findIndex((call) => call.id === callToDelete.id);
          if (runIdx >= 0) {
            this.onRunsChanged.value.splice(runIdx, 1);
            this.onRunsChanged.next(this.onRunsChanged.value);
          }
          this.onDelete.next(callToDelete);
        });
      });

    const anyMetadataEditSub = merge(...this.cards.map((card) => card.onMetadataEdit))
      .subscribe((editedCall) => {
        this.onMetadataEdit.next(editedCall);
      });

    const favoritesChanged = merge(...this.cards.map((card) => card.onFavoriteChanged))
      .subscribe(() => this.onFilteringChanged.next(this.onFilteringChanged.value));

    this.subs.push(anyClickSub, anySelectSub, anyUnselectSub, anyDeleteSub, anyMetadataEditSub, favoritesChanged);
  }

  private redrawSelectionState() {
    const currentSelectedSet = this._onSelectedChanged.value;

    const newTotal = ui.span([`Selected: ${currentSelectedSet.size}`], {style: {'align-self': 'center'}});
    this.total.replaceWith(newTotal);
    this.total = newTotal;

    if (currentSelectedSet.size < 2)
      this.compareIcon.classList.add('hp-disabled');
    else
      this.compareIcon.classList.remove('hp-disabled');

    if (currentSelectedSet.size === 0)
      this.trashIcon.classList.add('hp-disabled');
    else
      this.trashIcon.classList.remove('hp-disabled');

    if (currentSelectedSet.size === 0) {
      ui.setDisplay(this.selectAllIcon, true);
      ui.setDisplay(this.unselectOnPartialIcon, false);
      ui.setDisplay(this.unselectAllIcon, false);
    } else {
      if (currentSelectedSet.size === this.onFilteredRunsChanged.value.length) {
        ui.setDisplay(this.selectAllIcon, false);
        ui.setDisplay(this.unselectOnPartialIcon, false);
        ui.setDisplay(this.unselectAllIcon, true);
      } else {
        ui.setDisplay(this.selectAllIcon, false);
        ui.setDisplay(this.unselectOnPartialIcon, true);
        ui.setDisplay(this.unselectAllIcon, false);
      }
    }
  }
}
