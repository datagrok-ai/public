import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import {Subject, BehaviorSubject, merge, fromEvent} from 'rxjs';
import {historyUtils} from '../../history-utils';
import {properUpdateIndicator} from '../../function-views/src/shared/utils';
import {HistoricalRunsDelete, HistoricalRunEdit} from './history-dialogs';
import {EXPERIMENTAL_TAG} from '../../shared-utils/consts';

class HistoricalRunCard extends DG.Widget {
  private _onClicked = new Subject<DG.FuncCall>();
  private _onSelect = new Subject<DG.FuncCall>();
  private _onUnselect = new Subject<DG.FuncCall>();
  private _onMetadataEdit = new Subject<DG.FuncCall>();
  private _onDelete = new Subject<DG.FuncCall>();

  private onFuncCallChanged = new BehaviorSubject<DG.FuncCall>(this.initialFuncCall);

  private addToSelected = ui.iconFA('square', (ev) => {
    ev.stopPropagation();
    this._onSelect.next(this.funcCall);
    this.setSelected();
  }, 'Select this run');

  private removeFromSelected = ui.iconFA('check-square', (ev) => {
    ev.stopPropagation();
    this._onUnselect.next(this.funcCall);
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

  public onClicked = this._onClicked.asObservable();
  public onSelect = this._onSelect.asObservable();
  public onUnselect = this._onUnselect.asObservable();
  public onMetadataEdit = this._onMetadataEdit.asObservable();
  public onDelete = this._onDelete.asObservable();

  public get funcCall() {
    return this.onFuncCallChanged.value;
  }

  constructor(
    private initialFuncCall: DG.FuncCall,
    private options?: {
      showEdit?: boolean,
      showDelete?: boolean,
      showFavorite?: boolean,
      showAuthorIcon?: boolean,
    },
  ) {
    super(ui.div());

    const onChangedSub = merge(this.onFuncCallChanged, this.onMetadataEdit)
      .subscribe((funccall) => this.refresh(funccall));

    const clickSub = fromEvent(this.root, 'click').subscribe(() => this._onClicked.next(initialFuncCall));
    this.subs.push(clickSub, onChangedSub);
  }

  async setFavorited() {
    properUpdateIndicator(this.root, true);
    return historyUtils.loadRun(this.funcCall.id, false).then((fullCall) => {
      fullCall.options['isFavorite'] = true;
      return historyUtils.saveRun(fullCall);
    }).then((fullCall) => {
      ui.setDisplay(this.addToFavorites, false);
      ui.setDisplay(this.unfavoriteIcon, true);

      this._onMetadataEdit.next(fullCall);
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
      ui.setDisplay(this.unfavoriteIcon, false);
      ui.setDisplay(this.addToFavorites, true);

      this._onMetadataEdit.next(fullCall);
    }).catch((err) => {
      grok.shell.error(err);
    }).finally(() => {
      properUpdateIndicator(this.root, false);
    });
  }

  setSelected() {
    ui.setDisplay(this.addToSelected, false);
    ui.setDisplay(this.removeFromSelected, true);
  }

  setUnselected() {
    ui.setDisplay(this.removeFromSelected, false);
    ui.setDisplay(this.addToSelected, true);
  }

  updateCard(updatedFunccall: DG.FuncCall) {
    this.onFuncCallChanged.next(updatedFunccall);
    this._onMetadataEdit.next(updatedFunccall);
  }

  async deleteCard(showDialog = true) {
    if (showDialog) {
      const deleteDialog = new HistoricalRunsDelete(new Set([this.funcCall]));

      const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async (setToDelete) => {
        properUpdateIndicator(this.root, true);
        Promise.all(
          wu(setToDelete.values()).map(async () => {
            return historyUtils.deleteRun(this.funcCall).then(() => {
              this._onDelete.next(this.funcCall);
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

        this._onDelete.next(this.funcCall);
      });
    }

    return this._onDelete.toPromise();
  }

  private async editCard() {
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
          this._onMetadataEdit.next(fullCall);

          onEditSub.unsubscribe();
        })
        .catch((err) => {
          grok.shell.error(err);
        }).finally(() => {
          properUpdateIndicator(this.root, false);
        });
    });
    editDialog.show({center: true, width: 500});

    return this._onMetadataEdit.toPromise();
  }

  private refresh(funcCall: DG.FuncCall) {
    const options = this.options;

    let authorIcon = null;
    if (this.options?.showAuthorIcon) {
      authorIcon = funcCall.author.picture as HTMLElement;
      authorIcon.style.width = '25px';
      authorIcon.style.height = '25px';
      authorIcon.style.fontSize = '20px';

      ui.bind(funcCall.author, authorIcon);
    }

    const experimentalTag = ui.iconFA('flask', null, 'Experimental run');
    experimentalTag.classList.add('fad', 'fa-sm');
    experimentalTag.classList.remove('fal');
    experimentalTag.style.marginLeft = '3px';
    const immutableTags = funcCall.options['immutable_tags'] as string[] | undefined;
    const cardLabel = ui.span([
      ui.label(
        funcCall.options['title'] ??
        funcCall.author?.friendlyName ??
        grok.shell.user.friendlyName, {style: {'color': 'var(--blue-1)'}},
      ),
      ...(immutableTags && immutableTags.includes(EXPERIMENTAL_TAG)) ?
        [experimentalTag]:[],
    ]);

    const editIcon = ui.iconFA('edit', (ev) => {
      ev.stopPropagation();
      this.editCard();
    }, 'Edit run metadata');
    editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
      ev.stopPropagation();
      this.deleteCard();
    }, 'Delete run');
    deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

    this.addToFavorites.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.unfavoriteIcon.classList.add('fas', 'hp-funccall-card-icon');

    if (funcCall.options['isFavorite']) {
      ui.setDisplay(this.addToFavorites, false);
      ui.setDisplay(this.unfavoriteIcon, true);
    } else {
      ui.setDisplay(this.unfavoriteIcon, false);
      ui.setDisplay(this.addToFavorites, true);
    }

    this.addToSelected.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');
    this.removeFromSelected.classList.add('hp-funccall-card-icon');
    ui.setDisplay(this.removeFromSelected, false);

    const dateStarted = new Date(funcCall.started.toString())
      .toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'});

    const card = ui.divH([
      ui.divH([
        ...(options?.showAuthorIcon) ? [authorIcon!]: [],
        ui.divV([
          cardLabel,
          ui.span([dateStarted]),
          ...(funcCall.options['description']) ? [ui.divText(funcCall.options['description'], 'description')]: [],
          ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
            [ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')))]:[],
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
      ...(funcCall.options['title']) ? {'Title': funcCall.options['title']}:{},
      ...(funcCall.options['description']) ? {'Description': funcCall.options['description']}:{},
      ...(funcCall.options['tags'] && funcCall.options['tags'].length > 0) ?
        {'Tags': ui.div(funcCall.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor')}:{},
      ...(funcCall.options['immutable_tags'] && funcCall.options['immutable_tags'].length > 0) ?
        {'Extra': ui.div(funcCall.options['immutable_tags']
          .map((tag: string) => ui.span([tag], 'd4-tag')), 'd4-tag-editor',
        )}:{},
    }));

    ui.empty(this.root);
    this.root.appendChild(card);
  }
}

export class HistoricalRunsList extends DG.Widget {
  public onComparisonCalled = new Subject<string[]>();

  public onMetadataEdit = new Subject<DG.FuncCall>();
  public onDelete = new Subject<DG.FuncCall>();
  public onClicked = new Subject<DG.FuncCall>();

  private cards = this.initialFuncCalls.map((funcCall) => new HistoricalRunCard(funcCall, this.options));

  private onSelectedRunsChanged = new BehaviorSubject<Set<DG.FuncCall>>(new Set());
  private onRunsChanged = new BehaviorSubject<DG.FuncCall[]>(this.initialFuncCalls);

  private actions = ui.box() as HTMLElement;

  public get selected() {
    return this.onSelectedRunsChanged.value;
  }

  public get onSelectedChanged() {
    return this.onSelectedRunsChanged.asObservable();
  }

  constructor(private initialFuncCalls: DG.FuncCall[], private options: {
      fallbackText?: string,
      showEdit?: boolean,
      showDelete?: boolean,
      showFavorite?: boolean,
      showAuthorIcon?: boolean,
      showCompare?: boolean,
    } = {showEdit: true, showDelete: true, showFavorite: true, showAuthorIcon: true, showCompare: true}) {
    super(ui.div([], {style: {'height': '100%', 'overflow-y': 'hidden'}}));

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

  async deleteItem(id: string) {
    const corrCards = this.cards.splice(this.cards.findIndex((card) => card.funcCall.id === id)!, 1);
    // We delete single item
    await corrCards[0].deleteCard(false);
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
            await this.deleteItem(funcCall.id);
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
    const t = ui.iconFA('square',
      () => this.onSelectedRunsChanged.next(new Set(this.onRunsChanged.value)),
      'Select all',
    );
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
    this.cards = newFuncCalls.map((funcCall) => new HistoricalRunCard(funcCall, this.options));

    if (newFuncCalls.length > 0) {
      this.root.appendChild(ui.divV([
        this.actions,
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

        this.onRunsChanged.value.splice(this.onRunsChanged.value.findIndex((call) => call.id === callToDelete.id), 1);
        this.onRunsChanged.next( this.onRunsChanged.value);

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
        ...this.options.showCompare ? [this.getCompareIcon]: [],
        this.getTrashIcon,
        currentSelectedSet.size === 0 ? this.selectAllIcon: this.unselectAllIcon,
      ]),
    ], {style: {
      'justify-content': 'space-between',
      'padding': '0px 10px 0px 9px',
    }}) as HTMLElement;

    return actionsSection;
  }
}
