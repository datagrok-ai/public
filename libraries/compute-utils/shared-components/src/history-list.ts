import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import $ from 'cash-dom';
import {Subject, BehaviorSubject} from 'rxjs';
import {historyUtils} from '../../history-utils';
import {HistoricalRunsDelete, HistoricalRunEdit} from './history-dialogs';
import {EXPERIMENTAL_TAG, storageName} from '../../shared-utils/consts';
import {ID_COLUMN_NAME} from './history-input';

const EXP_COLUMN_NAME = 'Source';
const FAVORITE_COLUMN_NAME = 'Is favorite';
const ACTIONS_COLUMN_NAME = 'Delete';

const AUTHOR_COLUMN_NAME = 'Author';
const STARTED_COLUMN_NAME = 'Started';
const TITLE_COLUMN_NAME = 'Title';
const DESC_COLUMN_NAME = 'Desc.';
const TAGS_COLUMN_NAME = 'Tags';

const getColumnName = (key: string) => {
  return `${key[0].toUpperCase()}${key.substring(1)}`;
};

const runCache = new DG.LruCache<string, DG.FuncCall>();

export class HistoricalRunsList extends DG.Widget {
  private storageName(runs: Map<string, DG.FuncCall>) {
    return `${storageName}_${[...runs.values()][0].func.name}_Fav`;
  }

  private onRunsChanged = new BehaviorSubject<Map<string, DG.FuncCall>>(this.initialRuns.reduce((acc, run) => {
    acc.set(run.id, run);
    return acc;
  }, new Map<string, DG.FuncCall>));
  private onRunsDfChanged = new BehaviorSubject<DG.DataFrame>(DG.DataFrame.fromColumns([
    DG.Column.bool(EXP_COLUMN_NAME, 0),
    ...this.options.isHistory ? [DG.Column.bool(FAVORITE_COLUMN_NAME, 0)]: [],
    ...this.options.showActions ? [DG.Column.string(ACTIONS_COLUMN_NAME, 0)]: [],
    DG.Column.dateTime(STARTED_COLUMN_NAME, 0),
    DG.Column.string(AUTHOR_COLUMN_NAME, 0),
    DG.Column.string(TAGS_COLUMN_NAME, 0),
    DG.Column.string(TITLE_COLUMN_NAME, 0),
    DG.Column.string(DESC_COLUMN_NAME, 0),
    ...this.visibleProps
      .map((column) => DG.Column.fromStrings(getColumnName(column), [])),
    DG.Column.fromStrings(ID_COLUMN_NAME, []),
  ]));

  private _historyFilters = DG.Viewer.filters(this.onRunsDfChanged.value, {title: 'Filters'});
  private _historyGrid = DG.Viewer.grid(
    this.onRunsDfChanged.value,
    {showRowHeader: false, showColumnGridlines: false, allowEdit: false},
  );

  public onComparisonCalled = new Subject<string[]>();
  public onMetadataEdit = new Subject<DG.FuncCall>();
  public onClicked = new Subject<DG.FuncCall>();
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

  private showEditDialog(funcCall: DG.FuncCall, isFavorite: boolean) {
    const editDialog = new HistoricalRunEdit(funcCall, isFavorite);

    const onEditSub = editDialog.onMetadataEdit.subscribe(async (editOptions) => {
      if (!this.options.isHistory)
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

  constructor(
    private readonly initialRuns: DG.FuncCall[],
    private visibleProps: string[] = [],
    private readonly options: {
      fallbackText?: string,
      showActions?: boolean,
      showBatchActions?: boolean,
      // Used in the only place to avoid bug
      isHistory?: boolean
    } = {
      showActions: true,
      showBatchActions: true,
      isHistory: false,
    }) {
    super(ui.div([], {style: {height: '100%', width: '100%'}}));

    this.root.appendChild(this.defaultText);
    $(this.showFiltersIcon.firstChild).addClass('fad');

    const batchActions = ui.divH([
      ui.divH([
        this.showFiltersIcon, this.showMetadataIcon.root, this.showInputsIcon.root,
      ], {style: {'gap': '5px', 'padding': '0px 5px'}}),
      ui.divH([
        this.trashIcon, this.compareIcon,
      ], {style: {'gap': '5px', 'padding': '0px 5px'}}),
    ], {style: {'justify-content': 'space-between'}});
    batchActions.style.setProperty('overflow-y', 'hidden', 'important');

    const gridWithControls = ui.divV([
      ...this.options.showBatchActions ? [batchActions]:[],
      this._historyGrid.root,
    ], 'ui-box');
    $(gridWithControls).removeClass('ui-div');
    const splitH = ui.splitH([
      gridWithControls,
      this._historyFilters.root,
    ], {style: {'height': '100%'}}, true);
    this.root.appendChild(splitH);
    ui.tools.waitForElementInDom(splitH).then((x) => {
      $(x).css('width', '100%');
      $(this._historyFilters.root).css('width', '300px');
    });

    ui.setDisplayAll([this._historyGrid.root, this._historyFilters.root,
      this.showMetadataIcon.root, this.showInputsIcon.root, this.trashIcon,
      this.compareIcon, this.showFiltersIcon], false);

    const listChangedSub = this.onRunsChanged.subscribe(async (runs) => {
      const newRuns = [...runs.values()];

      const extractStringValue = (run: DG.FuncCall, key: string) => {
        if (key === AUTHOR_COLUMN_NAME) return run.author?.friendlyName ?? grok.shell.user.friendlyName;

        const val =
        (run as any)[key] ??
        run.inputs[key] ??
        run.outputs[key] ??
        run.options[key] ??
        null;

        return val?.toString() ?? '';
      };

      const getColumnByProp = (prop: DG.Property) => {
        if (prop.propertyType === DG.TYPE.DATE_TIME) {
          return DG.Column.dateTime(prop.caption ?? getColumnName(prop.name), newRuns.length)
            .init((idx) => (<any>window).grok_DayJs_To_DateTime(newRuns[idx].inputs[prop.name]));
        }

        return DG.Column.fromType(
          prop.propertyType as any,
          prop.caption ?? getColumnName(prop.name),
          newRuns.length,
        ).init((idx) => extractStringValue(newRuns[idx], prop.name));
      };

      const getColumnByName = (key: string) => {
        if (key === STARTED_COLUMN_NAME) {
          return DG.Column.dateTime(getColumnName(key), newRuns.length)
            .init((idx) => (<any>window).grok_DayJs_To_DateTime(newRuns[idx].started));
        }

        return DG.Column.fromStrings(getColumnName(key), newRuns.map((run) => extractStringValue(run, key)));
      };

      if (runs.size > 0) {
        const favoritesRecord: Record<string, string> =
          await grok.dapi.userDataStorage.get(this.storageName(runs)) ?? {};
        const favorites = Object.keys(favoritesRecord);

        const func = [...this.runs.values()][0].func;

        if (this.visibleProps.length === 0) {
          this.visibleProps = func.inputs
            .filter((input) => input.propertyType !== DG.TYPE.DATA_FRAME)
            .map((prop) => prop.name);
        }

        const newRunsGridDf = DG.DataFrame.fromColumns([
          DG.Column.string(EXP_COLUMN_NAME, newRuns.length).init((idx) => {
            const immutableTags = newRuns[idx].options['immutable_tags'] as string[];
            return immutableTags && immutableTags.includes(EXPERIMENTAL_TAG) ? 'Experimental': 'Simulated';
          }),
          ...this.options.isHistory ?
            [
              DG.Column.bool(FAVORITE_COLUMN_NAME, newRuns.length)
                .init((idx) => favorites.includes(newRuns[idx].id)),
            ]: [],
          ...this.options.showActions ? [DG.Column.string(ACTIONS_COLUMN_NAME, newRuns.length).init('')]: [],
          getColumnByName(STARTED_COLUMN_NAME),
          getColumnByName(AUTHOR_COLUMN_NAME),
          DG.Column.string(TAGS_COLUMN_NAME, newRuns.length).init((idx) =>
            newRuns[idx].options['tags'] ? newRuns[idx].options['tags'].join(','): '',
          ),
          DG.Column.string(TITLE_COLUMN_NAME, newRuns.length).init((idx) => newRuns[idx].options['title']),
          DG.Column.string(DESC_COLUMN_NAME, newRuns.length).init((idx) => newRuns[idx].options['description']),
          ...this.visibleProps.map((key) => getColumnByProp(func.inputs.find((prop) => prop.name === key)!)),
          DG.Column.fromStrings(ID_COLUMN_NAME, newRuns.map((newRun) => newRun.id)),
        ]);

        $(this.defaultText).hide();
        ui.setDisplayAll([this._historyGrid.root, this._historyFilters.root,
          this.showInputsIcon.root, this.showMetadataIcon.root, this.trashIcon,
          this.compareIcon, this.showFiltersIcon], true);

        this.onRunsDfChanged.next(newRunsGridDf);
      } else {
        $(this.defaultText).show();
        ui.setDisplayAll([this._historyGrid.root, this._historyFilters.root,
          this.showInputsIcon.root, this.showMetadataIcon.root, this.trashIcon,
          this.compareIcon, this.showFiltersIcon], false);
      }
    });

    this.onRunsDfChanged.subscribe((newRunsGridDf) => {
      this._historyFilters.dataFrame = newRunsGridDf;
      this._historyGrid.dataFrame = newRunsGridDf;

      this._historyGrid.dataFrame.getCol(TAGS_COLUMN_NAME).setTag(DG.TAGS.MULTI_VALUE_SEPARATOR, ',');

      this.styleHistoryGrid();
      this.styleHistoryFilters();

      newRunsGridDf.onCurrentRowChanged.subscribe(() => {
        this._historyGrid.dataFrame.rows.select(() => false);
        // Filtering out header clicks
        if (this._historyGrid.dataFrame.currentRowIdx >= 0)
          this.onClicked.next(this.runs.get(this._historyGrid.dataFrame.currentRow.get(ID_COLUMN_NAME)));
      });

      newRunsGridDf.onSelectionChanged.subscribe(() => {
        if (this._historyGrid.dataFrame.selection.trueCount > 0)
          this._historyGrid.setOptions({'showCurrentRowIndicator': false});
        else
          this._historyGrid.setOptions({'showCurrentRowIndicator': true});

        const selectedCalls: DG.FuncCall[] = [...newRunsGridDf.selection
          .getSelectedIndexes()]
          .map((idx) => {
            const call = this.runs.get(newRunsGridDf.rows.get(idx).get(ID_COLUMN_NAME))!;
            return call;
          });
        this._onSelectedChanged.next(new Set(selectedCalls));

        this.redrawSelectionState();
      });
    });

    this.subs.push(listChangedSub);
  }

  async deleteRun(id: string) {
    return historyUtils.loadRun(id, true)
      .then((loadedRun) => {
        return [
          (this.options.isHistory) ? historyUtils.deleteRun(loadedRun): Promise.resolve(),
          loadedRun,
        ] as const;
      })
      .then(([, loadedRun]) => {
        this.onRunsChanged.value.delete(id);
        this.onRunsChanged.next(this.onRunsChanged.value);
        if (this.options.isHistory) this.onDelete.next(loadedRun);
      })
      .catch((e) => {
        grok.shell.error(e);
      });
  }

  editRun(editedCall: DG.FuncCall) {
    this.showEditDialog(editedCall, false);
  }

  updateRun(updatedRun: DG.FuncCall) {
    this.onRunsChanged.value.set(updatedRun.id, updatedRun);
    this.onRunsChanged.next(this.onRunsChanged.value);
  }

  addRun(newRun: DG.FuncCall) {
    this.onRunsChanged.value.set(newRun.id, newRun);
    this.onRunsChanged.next(this.onRunsChanged.value);
  }

  updateRuns(newRuns: DG.FuncCall[]) {
    this.onRunsChanged.next(newRuns.reduce((acc, run) => {
      acc.set(run.id, run);
      return acc;
    }, new Map<string, DG.FuncCall>));
  }

  private trashIcon = ui.div(ui.iconFA('trash-alt', () => {
    const selection = this._historyGrid.dataFrame.selection;
    const setToDelete = new Set<DG.FuncCall>();
    for (let i = 0; i < selection.length; i++) {
      const bit = this._historyGrid.dataFrame.selection.get(i);
      if (bit)
        setToDelete.add(this.runs.get(this._historyGrid.dataFrame.getCol(ID_COLUMN_NAME).get(i))!);
    }

    const deleteDialog = new HistoricalRunsDelete(setToDelete);

    const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async (setToDelete) => {
      await Promise.all(
        wu(setToDelete.values()).map(async (funcCall) => {
          await this.deleteRun(funcCall.id);

          return Promise.resolve();
        }));

      onDeleteSub.unsubscribe();
    });

    deleteDialog.show({center: true, width: 500, modal: true});
  }, 'Delete selected runs'), {style: {'padding-top': '5px'}});

  private compareIcon = ui.div(ui.iconFA('exchange', async () => {
    const selection = this._historyGrid.dataFrame.selection;
    const setToCompare = [];
    for (let i = 0; i < selection.length; i++) {
      const bit = this._historyGrid.dataFrame.selection.get(i);
      if (bit)
        setToCompare.push(this._historyGrid.dataFrame.getCol(ID_COLUMN_NAME).get(i));
    }

    this.onComparisonCalled.next(setToCompare);
  }, 'Compare selected runs'), {style: {'padding-top': '5px'}});


  private isFilterHidden = false;
  private showFiltersIcon = ui.div(ui.iconFA('filter', () => {
    if (this.isFilterHidden) {
      $(this.showFiltersIcon.firstChild).addClass('fad');
      $(this.showFiltersIcon.firstChild).removeClass('fal');
      $(this._historyFilters.root.parentElement).css('width', '100%');
      $(this._historyFilters.root).css('width', '300px');
      ui.setDisplay(this._historyFilters.root, true);
      this.isFilterHidden = false;
    } else {
      $(this.showFiltersIcon.firstChild).addClass('fal');
      $(this.showFiltersIcon.firstChild).removeClass('fad');
      ui.setDisplay(this._historyFilters.root, false);
      this.isFilterHidden = true;
    }
  }, 'Toggle filters'), {style: {'padding-top': '4px'}});

  private showMetadataIcon = ui.switchInput('Metadata', true, (newValue: boolean) => {
    this.styleHistoryGrid();
    this.styleHistoryFilters();
  });

  private showInputsIcon = ui.switchInput('Inputs', false, async (newValue: boolean) => {
    if (this.runs.size === 0) return;

    if (newValue && this.options.isHistory) {
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
  });

  private defaultText = ui.divV([
    ui.element('div', 'splitbar-horizontal'),
    ui.divText(this.options?.fallbackText ?? 'No historical runs found', 'hp-no-elements-label'),
  ]);

  private styleHistoryGrid() {
    if (this.runs.size > 0) {
      const func = [...this.runs.values()][0].func;

      const actionsCol = this._historyGrid.columns.byName(ACTIONS_COLUMN_NAME);
      if (actionsCol) {
        actionsCol.cellType = 'html';
        actionsCol.width = 35;
      }

      const favCol = this._historyGrid.columns.byName(FAVORITE_COLUMN_NAME);
      if (favCol) {
        favCol.cellType = 'html';
        favCol.width = 20;
      }
      this._historyGrid.columns.byName(EXP_COLUMN_NAME)!.cellType = 'html';
      this._historyGrid.columns.byName(EXP_COLUMN_NAME)!.width = 20;
      this._historyGrid.columns.byName(TAGS_COLUMN_NAME)!.cellType = 'html';
      this._historyGrid.columns.byName(TAGS_COLUMN_NAME)!.width = 90;
      this._historyGrid.columns.byName(STARTED_COLUMN_NAME)!.width = 110;
      this._historyGrid.setOptions({
        'showCurrentCellOutline': false,
        'allowEdit': false,
        'allowBlockSelection': false,
      });
      for (let i = 0; i < this._historyGrid.columns.length; i++) {
        const col = this._historyGrid.columns.byIndex(i)?.column;
        if (col && col.type === DG.TYPE.DATE_TIME)
          this._historyGrid.columns.byIndex(i)!.format = 'MMM d HH:mm';
      }

      this._historyGrid.onCellPrepare((cell) => {
        if (cell.tableColumn?.name === ACTIONS_COLUMN_NAME) {
          if (cell.isColHeader) {
            cell.customText = '';
            return;
          }

          const run = this.runs.get(cell.tableRow?.get(ID_COLUMN_NAME))!;

          cell.element = ui.divH([
            ui.iconFA('trash', () => {
              const setToDelete= new Set([this.runs.get(cell.tableRow!.get(ID_COLUMN_NAME))!]);
              const deleteDialog = new HistoricalRunsDelete(setToDelete);

              const onDeleteSub = deleteDialog.onFuncCallDelete.subscribe(async () => {
                await Promise.all(
                  wu(setToDelete.values()).map(async (funcCall) => {
                    await this.deleteRun(funcCall.id);

                    return Promise.resolve();
                  }));
                onDeleteSub.unsubscribe();
              });
              deleteDialog.show({center: true, width: 500});
            }, 'Remove run from history'),
            ui.iconFA(
              'edit',
              () => this.showEditDialog(
                run,
                cell.tableRow!.get(FAVORITE_COLUMN_NAME)!,
              ),
              'Edit run metadata',
            ),
          ], {style: {'padding': '6px 0px', 'gap': '6px', 'justify-content': 'space-between'}});
        }

        if (cell.tableColumn?.name === FAVORITE_COLUMN_NAME) {
          if (cell.isColHeader) {
            cell.customText = '';
            return;
          }

          const run = this.runs.get(cell.tableRow?.get(ID_COLUMN_NAME))!;

          const unfavoriteIcon =
            ui.iconFA('star', () => this.saveIsFavorite(run, false).then(() => this.updateRun(run)), 'Unfavorite');
          unfavoriteIcon.classList.add('fas');

          cell.element = ui.div(
            cell.cell.value ?
              unfavoriteIcon:
              ui.iconFA('star', () => this.saveIsFavorite(run, true).then(() => this.updateRun(run)), 'Favorite'),
            {style: {'padding': '5px 0px'}});
        }

        if (cell.tableColumn?.name === EXP_COLUMN_NAME) {
          if (cell.isColHeader)
            cell.customText = '';

          const experimentalTag = ui.iconFA('flask', null, 'Experimental run');
          experimentalTag.classList.add('fad', 'fa-sm');
          experimentalTag.classList.remove('fal');

          cell.element = cell.cell.value && cell.cell.value === 'Experimental' ?
            ui.div(experimentalTag, {style: {'padding': '5px'}}): ui.div();
        }

        if (cell.tableColumn?.name === TAGS_COLUMN_NAME) {
          if (cell.isColHeader)
            return;

          const tags = cell.cell.value.length > 0 ? ui.div((cell.cell.value as string | null)?.split(',').map(
            (tag: string) => ui.span([tag], 'd4-tag')),
          'd4-tag-editor'): ui.div();
          $(tags).css({'padding': '3px', 'background-color': 'transparent'});
          cell.element = tags;
        }
      });

      const tagCol = this._historyGrid.dataFrame.getCol(TAGS_COLUMN_NAME);
      this._historyGrid.columns.setVisible([
        ...this.showMetadataIcon.value ? [EXP_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [FAVORITE_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [ACTIONS_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [STARTED_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [AUTHOR_COLUMN_NAME]: [],
        ...tagCol.stats.missingValueCount < tagCol.length && this.showMetadataIcon.value ? [TAGS_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [TITLE_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [DESC_COLUMN_NAME]: [],
        ...this.showInputsIcon.value ? this.visibleProps
          .map((visibleCol) => func.inputs.find((input) => input.name === visibleCol)!)
          .filter((prop) => {
            return !this.options.isHistory ||
              this._historyGrid.dataFrame.getCol(prop.caption ?? getColumnName(prop.name)).categories.length > 1;
          })
          .map((prop) => prop.caption ?? getColumnName(prop.name)): [],
      ]);

      this._historyGrid.root.style.height = '100%';
    }
  }

  private styleHistoryFilters() {
    if (this.runs.size > 0) {
      const func = [...this.runs.values()][0].func;

      const tagCol = this._historyGrid.dataFrame.getCol(TAGS_COLUMN_NAME);
      this._historyFilters.setOptions({columnNames: [
        ...this.showMetadataIcon.value &&
        this._historyGrid.dataFrame.getCol(EXP_COLUMN_NAME).categories.length > 1 ? [EXP_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value &&
        (this._historyGrid.dataFrame.col(FAVORITE_COLUMN_NAME)?.categories.length ?? 0) > 1 ?
          [FAVORITE_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value ? [STARTED_COLUMN_NAME]:[],
        ...this.showMetadataIcon.value &&
        this._historyGrid.dataFrame.getCol(AUTHOR_COLUMN_NAME).categories.length > 1 ? [AUTHOR_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value &&
        tagCol.stats.missingValueCount < tagCol.length ? [TAGS_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value &&
        this._historyGrid.dataFrame.getCol(TITLE_COLUMN_NAME).categories.length > 1 ? [TITLE_COLUMN_NAME]: [],
        ...this.showMetadataIcon.value &&
        this._historyGrid.dataFrame.getCol(DESC_COLUMN_NAME).categories.length > 1 ? [DESC_COLUMN_NAME]: [],
        ...this.showInputsIcon.value ? this.visibleProps
          .map((visibleCol) => func.inputs.find((input) => input.name === visibleCol)!)
          .filter((prop) => {
            return !this.options.isHistory ||
            this._historyGrid.dataFrame.getCol(prop.caption ?? getColumnName(prop.name)).categories.length > 1;
          })
          .map((prop) => prop.caption ?? getColumnName(prop.name)): [],
      ]});
      $(this._historyFilters.root).find('.d4-filter-group-header').hide();

      this._historyFilters.root.style.height = '100%';
    }
  }

  private redrawSelectionState() {
    const currentSelectedCount = this._historyGrid.dataFrame.selection.trueCount;

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
