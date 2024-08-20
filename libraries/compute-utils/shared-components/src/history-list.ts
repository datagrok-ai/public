import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import $ from 'cash-dom';
import {Subject, BehaviorSubject} from 'rxjs';
import {take} from 'rxjs/operators';
import {historyUtils} from '../../history-utils';
import {HistoricalRunsDelete, HistoricalRunEdit} from './history-dialogs';
import {ACTIONS_COLUMN_NAME, AUTHOR_COLUMN_NAME,
  COMPLETE_COLUMN_NAME,
  DESC_COLUMN_NAME, EXPERIMENTAL_TAG, EXP_COLUMN_NAME,
  FAVORITE_COLUMN_NAME, HistoryOptions, STARTED_COLUMN_NAME, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME
  , storageName} from '../../shared-utils/consts';
import {ID_COLUMN_NAME} from './history-input';
import {camel2title, extractStringValue, getColumnName, getMainParams, getRunsDfFromList, getStartedOrNull, getVisibleProps, styleHistoryFilters, styleHistoryGrid} from '../../shared-utils/utils';
import {getStarted} from '../../function-views/src/shared/utils';
import dayjs from 'dayjs';
import utc from 'dayjs/plugin/utc';

dayjs.extend(utc);

const runCache = new DG.LruCache<string, DG.FuncCall>();

export class HistoricalRunsList extends DG.Widget {
  private onRunsChanged = new BehaviorSubject<Map<string, DG.FuncCall>>(this.initialRuns.reduce((acc, run) => {
    acc.set(run.id, run);
    return acc;
  }, new Map<string, DG.FuncCall>));
  private _onRunsDfChanged = new BehaviorSubject<DG.DataFrame>(
    DG.DataFrame.fromColumns([
      DG.Column.bool(EXP_COLUMN_NAME, 0),
      ...this.options?.isHistory ? [DG.Column.bool(FAVORITE_COLUMN_NAME, 0)]: [],
      ...this.options?.showActions ? [DG.Column.string(ACTIONS_COLUMN_NAME, 0)]: [],
      DG.Column.bool(COMPLETE_COLUMN_NAME, 0),
      DG.Column.dateTime(STARTED_COLUMN_NAME, 0),
      DG.Column.string(AUTHOR_COLUMN_NAME, 0),
      DG.Column.string(TAGS_COLUMN_NAME, 0),
      DG.Column.string(TITLE_COLUMN_NAME, 0),
      DG.Column.string(DESC_COLUMN_NAME, 0),
      DG.Column.fromStrings(ID_COLUMN_NAME, []),
    ]));

  public get onRunsDfChanged() {
    return this._onRunsDfChanged.asObservable();
  }

  private _historyFilters = DG.Viewer.filters(this._onRunsDfChanged.value,
    {title: 'Filters', columnNames: [] as string[]});
  private _historyGrid = DG.Viewer.grid(
    this._onRunsDfChanged.value,
    {showRowHeader: false, showColumnGridlines: false, allowEdit: false},
  );

  public onComparisonCalled = new Subject<string[]>();
  public onMetadataEdit = new Subject<DG.FuncCall>();
  private _onChosen = new BehaviorSubject<DG.FuncCall | null>(null);
  public get onChosen() {
    return this._onChosen.asObservable();
  }

  public get chosen() {
    return this._onChosen.value;
  }

  public set chosen(val: DG.FuncCall | null) {
    this.chosenById = val?.id ?? null;
  }

  public set chosenById(id: string | null) {
    if (!id) {
      this.currentDf.currentRowIdx = -1;
      return;
    }

    for (let i = 0; i < this.currentDf.rowCount; i++) {
      if (this.currentDf.getCol(ID_COLUMN_NAME).get(i) === id) {
        this.currentDf.currentRowIdx = i;
        break;
      }
    }
  }

  public onDelete = new Subject<DG.FuncCall>();

  private _onSelectedChanged = new BehaviorSubject<Set<DG.FuncCall>>(new Set([]));

  public get selected() {
    return this._onSelectedChanged.value;
  }

  public get runs() {
    return this.onRunsChanged.value;
  }

  public get onSelectedChanged() {
    return this._onSelectedChanged.asObservable();
  }

  private async saveIsFavorite(funcCall: DG.FuncCall, isFavorite: boolean) {
    const favStorageName = `${storageName}_${funcCall.func.name}_Fav`;

    if (isFavorite)
      return grok.dapi.userDataStorage.postValue(favStorageName, funcCall.id, '');
    else
      return grok.dapi.userDataStorage.remove(favStorageName, funcCall.id);
  }

  private getRunByIdx(idx: number) {
    return this.runs.get(this.currentDf.get(ID_COLUMN_NAME, idx));
  }

  private showEditDialog(funcCall: DG.FuncCall, isFavorite: boolean) {
    const editDialog = new HistoricalRunEdit(funcCall, isFavorite);

    const onEditSub = editDialog.onMetadataEdit.subscribe(async (editOptions) => {
      if (!this.options?.isHistory)
        this.updateRun(funcCall);
      else {
        return ((editOptions.favorite !== 'same') ?
          this.saveIsFavorite(funcCall, (editOptions.favorite === 'favorited')) :
          Promise.resolve())
          .then(() => historyUtils.loadRun(funcCall.id, false))
          .then((fullCall) => {
            if (editOptions.title) fullCall.options['title'] = editOptions.title;
            if (editOptions.description) fullCall.options['description'] = editOptions.description;
            if (editOptions.tags) fullCall.options['tags'] = editOptions.tags;

            return [historyUtils.saveRun(fullCall), fullCall] as const;
          })
          .then(([, fullCall]) => {
            this.updateRun(fullCall);

            this.onMetadataEdit.next(fullCall);

            onEditSub.unsubscribe();
          })
          .catch((err) => {
            grok.shell.error(err);
          });
      }
    });
    editDialog.show({center: true, width: 500});
  }

  private layout = ui.splitH([], null, true);

  private defaultGridText = ui.divText(
    this.options?.fallbackText ?? 'No historical runs found', 'hp-no-elements-label',
  );

  private defaultFiltersText = ui.divText('No filters to show', 'hp-no-elements-label');

  private filterWithText = ui.divH([
    this.defaultFiltersText,
    this._historyFilters.root,
  ], 'ui-box');

  private visibleProps(func: DG.Func) {
    return getVisibleProps(func, this.options)
  };

  constructor(
    private readonly initialRuns: DG.FuncCall[],
    private options?: HistoryOptions) {
    super(ui.div([], {style: {height: '100%', width: '100%'}}));

    const batchActions = ui.divH([
      ui.divH([
        ...options?.isHistory ? [this.toggleCompactMode]: [],
        this.showFiltersIcon, this.showMetadataIcon.root, this.showInputsIcon.root,
      ], {style: {'gap': '5px', 'padding': '5px'}}),
      ui.divH(
        this.options?.showBatchActions ? [this.trashIcon, this.compareIcon]: [],
        {style: {'gap': '5px', 'padding': '5px'}}),
    ], {style: {'justify-content': 'space-between', 'padding': '0px 5px'}});
    batchActions.style.setProperty('overflow-y', 'hidden', 'important');

    const gridWithControls = ui.divV([
      this.defaultGridText,
      batchActions,
      this._historyGrid.root,
    ], 'ui-box');
    $(gridWithControls).removeClass('ui-div');
    $(gridWithControls).css('width', '100%');

    this.root.appendChild(this.layout);

    ui.setDisplayAll([this._historyGrid.root, this._historyFilters.root,
      this.toggleCompactMode, this.trashIcon,
      this.compareIcon, this.showFiltersIcon], false);

    this.setGridCellRendering();

    const listChangedSub = this.onRunsChanged.subscribe(async (runs) => {
      if (runs.size > 0) {
        const newRunsDf = await getRunsDfFromList(
          runs, 
          this.options
        );
        const func = runs.values().next().value.func as DG.Func;

        ui.setDisplayAll([this.defaultGridText, this.defaultFiltersText], false);
        ui.setDisplayAll([
          this._historyGrid.root,
          this._historyFilters.root,
          this.toggleCompactMode,
          this.showMetadataIcon.root,
          this.trashIcon,
          this.compareIcon,
          this.showFiltersIcon,
        ], true);
        ui.setDisplay(this.showInputsIcon.root, getVisibleProps(func, this.options).length > 0);

        this._onRunsDfChanged.next(newRunsDf);
      } else {
        ui.setDisplayAll([this.defaultGridText, this.defaultFiltersText], true);
        ui.setDisplayAll([
          this._historyGrid.root,
          this._historyFilters.root,
          this.showMetadataIcon.root,
          this.showInputsIcon.root,
          this.toggleCompactMode,
          this.trashIcon,
          this.compareIcon,
          this.showFiltersIcon,
        ], false);
      }
    });

    const runDfChangedSub = this._onRunsDfChanged.subscribe((newRunsGridDf) => {
      this._historyFilters.dataFrame = newRunsGridDf;
      this._historyGrid.dataFrame = newRunsGridDf;

      this.currentDf.getCol(TAGS_COLUMN_NAME).meta.multiValueSeparator = ',';

      this.styleHistoryGrid();
      this.styleHistoryFilters();

      const currentDfSub = newRunsGridDf.onCurrentRowChanged.subscribe(() => {
        this.currentDf.rows.select(() => false);
        this._onChosen.next(this.getRunByIdx(this.currentDf.currentRowIdx) ?? null);
      });

      const selectionSub = newRunsGridDf.onSelectionChanged.subscribe(() => {
        const selectedCalls: DG.FuncCall[] = [...newRunsGridDf.selection
          .getSelectedIndexes()]
          .map((idx) => this.getRunByIdx(idx)!);
        this._onSelectedChanged.next(new Set(selectedCalls));

        this.redrawSelectionState();
      });

      this.redrawSelectionState();

      this.subs.push(currentDfSub, selectionSub);
    });

    const filterSub = this._isFilterHidden.subscribe((isHidden) => {
      if (!isHidden) {
        $(this.showFiltersIcon.firstChild).addClass('fad');
        $(this.showFiltersIcon.firstChild).removeClass('fal');
        $(this.filterWithText).css(this.compactMode ? 'height': 'width', '320px');
        $(this.filterWithText).css(this.compactMode ? 'width': 'height', '100%');
        ui.setDisplay(this.filterWithText, true);
      } else {
        $(this.showFiltersIcon.firstChild).addClass('fal');
        $(this.showFiltersIcon.firstChild).removeClass('fad');
        ui.setDisplay(this.filterWithText, false);
      }
    });

    const compactModeSub = this._compactMode.subscribe((newValue) => {
      const content = [
        gridWithControls,
        this.filterWithText,
      ];
      $(this.filterWithText).removeClass('ui-div');
      const styles = {style: {'width': '100%', 'height': '100%'}};

      $(this.toggleCompactMode.firstChild).addClass(newValue ? 'fa-expand-alt': 'fa-compress-alt');
      $(this.toggleCompactMode.firstChild).removeClass(newValue ? 'fa-compress-alt': 'fa-expand-alt');

      const newLayout = newValue ?
        ui.splitV(content, styles, true) :
        ui.splitH(content, styles, true);
      this.layout.replaceWith(newLayout);
      this.layout = newLayout;

      ui.tools.waitForElementInDom(this.layout).then((x) => {
        $(x).css(this.compactMode ? 'height': 'width', '100%');

        $(this._historyFilters.root).css({
          ...this.compactMode ? {'width': '100%'} : {'height': '100%'},
        });
        $(this.filterWithText).css({
          ...this.compactMode ? {'height': '310px', 'width': '100%'} : {'width': '310px', 'height': '100%'},
        });
        $(gridWithControls).css({
          ...this.compactMode ? {'width': '100%'} : {'height': '100%'},
        });

        setTimeout(() => {
          this.styleHistoryGrid();
          this.styleHistoryFilters();
        }, 100);
      });
    });

    this.subs.push(listChangedSub, compactModeSub, runDfChangedSub, filterSub);
  }

  private isFavoriteByIndex(idx: number) {
    return this.currentDf.get(FAVORITE_COLUMN_NAME, idx);
  }

  private setGridCellRendering() {
    this._historyGrid.onCellPrepare((cell) => {
      if (cell.isColHeader && cell.tableColumn?.name &&
          [ACTIONS_COLUMN_NAME, EXP_COLUMN_NAME, FAVORITE_COLUMN_NAME].includes(cell.tableColumn.name))
        cell.customText = '';

      if (cell.isColHeader)
        return;

      if (cell.tableColumn?.name === ACTIONS_COLUMN_NAME) {
        cell.customText = '';
        const run = this.getRunByIdx(cell.tableRowIndex!)!;

        cell.element = ui.divH([
          ui.iconFA('trash', () => {
            const setToDelete = new Set([this.getRunByIdx(cell.tableRowIndex!)!]);
            const deleteDialog = new HistoricalRunsDelete(setToDelete);

            deleteDialog.onFuncCallDelete.pipe(
              take(1),
            ).subscribe(async () => {
              ui.setUpdateIndicator(this.root, true);
              try {
                await Promise.all(
                  wu(setToDelete.values()).map(async (funcCall) => {
                    await this.deleteRun(funcCall.id);

                    return Promise.resolve();
                  }));
              } catch (e: any) {
                grok.shell.error(e);
              } finally {
                ui.setUpdateIndicator(this.root, false);
              }
            });
            deleteDialog.show({center: true, width: 500});
          }, 'Remove run from history'),
          ui.iconFA(
            'edit',
            () => this.showEditDialog(
              run,
              this.isFavoriteByIndex(cell.tableRowIndex!),
            ),
            'Edit run metadata',
          ),
        ], {style: {'padding': '6px 0px', 'gap': '6px', 'justify-content': 'space-between'}});
      }

      if (cell.tableColumn?.name === FAVORITE_COLUMN_NAME) {
        cell.customText = '';
        const run = this.getRunByIdx(cell.tableRowIndex!)!;

        const unfavoriteIcon =
          ui.iconFA('star', () => this.saveIsFavorite(run, false).then(() => this.updateRun(run)), 'Unfavorite');
        $(unfavoriteIcon).addClass('fas');

        cell.element = ui.div(
          cell.cell.value ?
            unfavoriteIcon :
            ui.iconFA('star', () => this.saveIsFavorite(run, true).then(() => this.updateRun(run)), 'Favorite'),
          {style: {'padding': '5px 0px'}});
      }

      if (cell.tableColumn?.name === EXP_COLUMN_NAME) {
        cell.customText = '';
        const experimentalTag = ui.iconFA('flask', null, 'Experimental run');
        $(experimentalTag).addClass('fad fa-sm');
        $(experimentalTag).removeClass('fal');

        cell.element = cell.cell.value && cell.cell.value === 'Experimental' ?
          ui.div(experimentalTag, {style: {'padding': '5px'}}) : ui.div();
      }

      if (cell.tableColumn?.name === TAGS_COLUMN_NAME) {
        cell.customText = '';
        const tags = cell.cell.value.length > 0 ? ui.div((cell.cell.value as string | null)?.split(',').map(
          (tag: string) => ui.span([tag], 'd4-tag')),
        'd4-tag-editor') : ui.div();
        $(tags).css({'padding': '3px', 'background-color': 'transparent'});
        cell.element = tags;
      }

      if (cell.tableColumn?.name === ID_COLUMN_NAME) {
        cell.customText = '';
        const run = this.getRunByIdx(cell.tableRowIndex!);

        if (!run) return;

        const authorIcon = run.author.picture as HTMLElement;
        $(authorIcon).css({'width': '25px', 'height': '25px', 'fontSize': '20px'});

        ui.bind(run.author, authorIcon);

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
          this.showEditDialog(run, this.isFavoriteByIndex(cell.tableRowIndex!));
        }, 'Edit run metadata');
        editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

        const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
          ev.stopPropagation();
          const setToDelete = new Set([this.getRunByIdx(cell.tableRowIndex!)!]);
          const deleteDialog = new HistoricalRunsDelete(setToDelete);
          deleteDialog.onFuncCallDelete.pipe(
            take(1),
          ) .subscribe(async () => {
            ui.setUpdateIndicator(this.root, true);
            try {
              await Promise.all(
                wu(setToDelete.values()).map(async (funcCall) => {
                  await this.deleteRun(funcCall.id);

                  return Promise.resolve();
                }));
            } catch (e: any) {
              grok.shell.error(e);
            } finally {
              ui.setUpdateIndicator(this.root, false);
            }
          });
          deleteDialog.show({center: true, width: 500});
        }, 'Delete run');
        deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');


        const unfavoriteIcon =
        ui.iconFA('star', (ev) => {
          ev.stopPropagation();
          this.saveIsFavorite(run, false).then(() => this.updateRun(run));
        }, 'Unfavorite');
        unfavoriteIcon.classList.add('fas', 'hp-funccall-card-icon');

        const addToFavorites = ui.iconFA('star',
          (ev) => {
            ev.stopPropagation();
            this.saveIsFavorite(run, true).then(() => this.updateRun(run));
          }, 'Favorite');
        addToFavorites.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

        if (this.isFavoriteByIndex(cell.tableRowIndex!)) {
          ui.setDisplay(addToFavorites, false);
          ui.setDisplay(unfavoriteIcon, true);
        } else {
          ui.setDisplay(unfavoriteIcon, false);
          ui.setDisplay(addToFavorites, true);
        }

        const dateStarted = getStarted(run);

        const card = ui.divH([
          ui.divH([
            authorIcon,
            ui.divV([
              cardLabel,
              ui.span([dateStarted]),
              ...(run.options['tags'] && run.options['tags'].length > 0) ?
                [ui.div(run.options['tags'].map((tag: string) => ui.span([tag], 'd4-tag')))]:[],
            ], 'hp-card-content'),
          ]),
          ui.divH([
            ...(this.options?.showActions) ? [unfavoriteIcon, addToFavorites, editIcon, deleteIcon]: [],
          ]),
        ], 'hp-funccall-card');

        const tableRowIndex = cell.tableRowIndex!;
        card.addEventListener('mouseover', () => {
          cell.grid.dataFrame.mouseOverRowIdx = tableRowIndex;
        });
        card.addEventListener('click', (e) => {
          if (e.shiftKey)
            cell.grid.dataFrame.selection.set(tableRowIndex, true);
          else if (e.ctrlKey)
            cell.grid.dataFrame.selection.set(tableRowIndex, false);
          else
            cell.grid.dataFrame.currentRowIdx = tableRowIndex;
        });
        cell.element = card;
      }
    });
  }

  async deleteRun(id: string) {
    return historyUtils.loadRun(id, true)
      .then(async (loadedRun) => {
        return [
          await (this.options?.isHistory ? historyUtils.deleteRun(loadedRun): Promise.resolve()),
          loadedRun,
        ] as const;
      })
      .then(([, loadedRun]) => {
        this.runs.delete(id);
        this.onRunsChanged.next(this.runs);
        if (this.options?.isHistory) this.onDelete.next(loadedRun);

        return loadedRun;
      })
      .then((loadedRun) => {
        this.saveIsFavorite(loadedRun, false);
      })
      .catch((e) => {
        grok.shell.error(e);
      });
  }

  editRun(editedCall: DG.FuncCall) {
    this.showEditDialog(editedCall, false);
  }

  updateRun(updatedRun: DG.FuncCall) {
    this.runs.set(updatedRun.id, updatedRun);
    this.onRunsChanged.next(this.onRunsChanged.value);
  }

  addRun(newRun: DG.FuncCall) {
    this.runs.set(newRun.id, newRun);
    this.onRunsChanged.next(this.onRunsChanged.value);
  }

  updateRuns(newRuns: DG.FuncCall[]) {
    this.onRunsChanged.next(newRuns.reduce((acc, run) => {
      acc.set(run.id, run);
      return acc;
    }, new Map<string, DG.FuncCall>));
  }

  private trashIcon = ui.div(ui.iconFA('trash-alt', () => {
    const selection = this.currentDf.selection.getSelectedIndexes();
    const setToDelete = new Set<DG.FuncCall>();
    selection.forEach((idx) => setToDelete.add(this.getRunByIdx(idx)!));

    const deleteDialog = new HistoricalRunsDelete(setToDelete);

    deleteDialog.onFuncCallDelete.pipe(
      take(1),
    ).subscribe(async (setToDelete) => {
      ui.setUpdateIndicator(this.root, true);
      try {
        await Promise.all(
          wu(setToDelete.values()).map(async (funcCall) => {
            await this.deleteRun(funcCall.id);

            return Promise.resolve();
          }));
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        ui.setUpdateIndicator(this.root, false);
      }
    });

    deleteDialog.show({center: true, width: 500, modal: true});
  }, 'Delete selected runs'), {style: {'padding-top': '5px'}});

  private compareIcon = ui.div(ui.iconFA('exchange', async () => {
    const selection = this.currentDf.selection.getSelectedIndexes();
    const setToCompare = [] as string[];
    selection.forEach((idx) => {
      setToCompare.push(this.currentDf.getCol(ID_COLUMN_NAME).get(idx));
    });

    this.onComparisonCalled.next(setToCompare);
  }, 'Compare selected runs'), {style: {'padding-top': '5px'}});


  private _isFilterHidden = new BehaviorSubject(false);
  public get isFilterHidden() {
    return this._isFilterHidden.value;
  }

  public set isFilterHidden(value: boolean) {
    this._isFilterHidden.next(value);
  }
  private showFiltersIcon = ui.div(ui.iconFA('filter', () => {
    this._isFilterHidden.next(!this._isFilterHidden.value);
  }, 'Toggle filters'), {style: {'padding-top': '4px'}});

  private _compactMode = new BehaviorSubject(this.options?.isHistory ?? false);
  public get compactMode() {
    return this._compactMode.value;
  }

  public set compactMode(value: boolean) {
    this._compactMode.next(value);
  }

  public get onCompactModeChanged() {
    return this._compactMode.asObservable();
  }

  private toggleCompactMode = ui.div(ui.iconFA('expand-alt', () => {
    this.compactMode = !this.compactMode;
  }, 'Toggle compact mode'), {style: {'padding-top': '4px'}});

  private showMetadataIcon = ui.input.toggle('Metadata', {value: true, onValueChanged: () => {
    this.styleHistoryGrid();
    this.styleHistoryFilters();
  }});

  private showInputsIcon = ui.input.toggle('Params', {value: false, onValueChanged: async () => {
    const newValue = this.showInputsIcon.value;
    if (this.runs.size === 0) return;

    if (newValue && this.options?.isHistory) {
      const fullCalls = await Promise.all(
        [...this.runs.values()].map(async (run) => {
          if (runCache.has(run.id))
            return runCache.get(run.id)!;
          else {
            const loadedRun = await historyUtils.loadRun(run.id, true);
            runCache.set(loadedRun.id, loadedRun);

            return loadedRun;
          }
        }));
      fullCalls.forEach((fullCall) => {
        this.runs.set(fullCall.id, fullCall);
      });
      this.onRunsChanged.next(this.runs);
    } else {
      this.styleHistoryGrid();
      this.styleHistoryFilters();
    }
  }});

  private styleHistoryGrid() {
    styleHistoryGrid(
      this._historyGrid,
      this._compactMode.value,
      this.showInputsIcon.value,
      this.showMetadataIcon.value,
      this.runs.values().next().value.func
    )
  }

  private get currentDf() {
    return this._onRunsDfChanged.value;
  }

  private styleHistoryFilters() {
    const func = this.runs.values().next().value.func as DG.Func | undefined;

    const hasValidFilters = styleHistoryFilters(
      this._historyFilters,
      this.showMetadataIcon.value,
      this.showInputsIcon.value,
      this.options?.isHistory ?? false,
      func
    )

    if (hasValidFilters) {
      $(this.filterWithText).css({paddingTop: '10px'});
      ui.setDisplay(this.defaultFiltersText, false);
    } else {
      ui.setDisplay(this.defaultFiltersText, true);
    }
  }

  private redrawSelectionState() {
    const currentSelectedCount = this.currentDf.selection.trueCount;

    if (currentSelectedCount < 2)
      (this.compareIcon.firstChild as HTMLElement).classList.add('hp-disabled');
    else
      (this.compareIcon.firstChild as HTMLElement).classList.remove('hp-disabled');

    if (currentSelectedCount === 0)
      (this.trashIcon.firstChild as HTMLElement).classList.add('hp-disabled');
    else
      (this.trashIcon.firstChild as HTMLElement).classList.remove('hp-disabled');
  }
}
