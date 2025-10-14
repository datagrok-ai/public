import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {ACTIONS_COLUMN_NAME, AUTHOR_COLUMN_NAME, COMPLETE_COLUMN_NAME, DESC_COLUMN_NAME, EXP_COLUMN_NAME, EXPERIMENTAL_TAG, FAVORITE_COLUMN_NAME, HistoryOptions, STARTED_COLUMN_NAME, HISTORY_SUPPORTED_COL_TYPES, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME, storageName, VERSION_COLUMN_NAME} from './consts';
import dayjs from 'dayjs';
import {ID_COLUMN_NAME} from './consts';
import {getStartedOrNull, isIncomplete, getStarted} from './utils';

export const getFavStorageName = (func: DG.Func) => {
  return `${storageName}_${func.name}_Fav`;
};

export const saveIsFavorite = async (funcCall: DG.FuncCall, isFavorite: boolean) => {
  const favStorageName = `${storageName}_${funcCall.func.name}_Fav`;

  if (isFavorite)
    return grok.dapi.userDataStorage.postValue(favStorageName, funcCall.id, '');
  else
    return grok.dapi.userDataStorage.remove(favStorageName, funcCall.id);
};

export const loadIsFavorite = async (funcCall: DG.FuncCall): Promise<boolean> => {
  const favStorageName = `${storageName}_${funcCall.func.name}_Fav`;
  const hasEntry = await grok.dapi.userDataStorage.getValue(favStorageName, funcCall.id, true);

  return (hasEntry === '') || false;
};

export const getMainParams = (func: DG.Func): string[] | null => {
  return func.options['mainParams'] ? JSON.parse(func.options['mainParams']): null;
};

export const getVisibleProps = (func: DG.Func, options?: HistoryOptions): string[] => {
  return options?.visibleProps ?? getMainParams(func) ?? func.inputs
    .filter((input) => HISTORY_SUPPORTED_COL_TYPES.includes(input.propertyType as any))
    .map((prop) => prop.name);
};

export const setGridCellRendering = (
  grid: DG.Grid,
  runs: Map<string, DG.FuncCall>,
  onEditClick: (cell: DG.GridCell) => void,
  onDeleteClick: (cell: DG.GridCell) => void,
  onFavoriteClick: (cell: DG.GridCell) => void,
  onUnfavoriteClick: (cell: DG.GridCell) => void,
  showActions: boolean,
) => {
  const getRunByIdx = (idx: number) => {
    return runs.get(grid.dataFrame.get(ID_COLUMN_NAME, idx));
  };

  const isFavoriteByIndex = (idx: number) => {
    return grid.dataFrame.get(FAVORITE_COLUMN_NAME, idx);
  };

  grid.onCellPrepare((cell) => {
    if (cell.isColHeader && cell.tableColumn?.name &&
        [ACTIONS_COLUMN_NAME, EXP_COLUMN_NAME, FAVORITE_COLUMN_NAME].includes(cell.tableColumn.name))
      cell.customText = '';

    if (cell.isColHeader)
      return;

    if (cell.tableColumn?.name === ACTIONS_COLUMN_NAME) {
      cell.customText = '';

      cell.element = ui.divH([
        ui.iconFA('trash', () => onDeleteClick(cell), 'Remove run from history'),
        ui.iconFA(
          'edit',
          () => onEditClick(cell),
          'Edit run metadata',
        ),
      ], {style: {padding: '6px', gap: '6px', justifyContent: 'space-between'}});
    }

    if (cell.tableColumn?.name === FAVORITE_COLUMN_NAME) {
      cell.customText = '';
      const unfavoriteIcon =
        ui.iconFA('star', () => onUnfavoriteClick(cell), 'Unfavorite');
      $(unfavoriteIcon).addClass('fas');

      cell.element = ui.div(
        cell.cell.value ?
          unfavoriteIcon :
          ui.iconFA('star', () => onFavoriteClick(cell), 'Favorite'),
        {style: {'padding': '6px'}});
    }

    if (cell.tableColumn?.name === EXP_COLUMN_NAME) {
      cell.customText = '';
      const experimentalTag = ui.iconFA('flask', null, 'Experimental run');
      $(experimentalTag).addClass('fad fa-sm');
      $(experimentalTag).removeClass('fal');

      cell.element = cell.cell.value && cell.cell.value === 'Experimental' ?
        ui.div(experimentalTag, {style: {'padding': '6px'}}) : ui.div();
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
      const run = getRunByIdx(cell.tableRowIndex!);

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
        onEditClick(cell);
      }, 'Edit run metadata');
      editIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

      const deleteIcon = ui.iconFA('trash-alt', async (ev) => {
        ev.stopPropagation();
        onDeleteClick(cell);
      }, 'Delete run');
      deleteIcon.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

      const unfavoriteIcon = ui.iconFA('star', (ev) => {
        ev.stopPropagation();
        onUnfavoriteClick(cell);
      }, 'Unfavorite');
      unfavoriteIcon.classList.add('fas', 'hp-funccall-card-icon');

      const addToFavorites = ui.iconFA('star',
        (ev) => {
          ev.stopPropagation();
          onFavoriteClick(cell);
        }, 'Favorite');
      addToFavorites.classList.add('hp-funccall-card-icon', 'hp-funccall-card-hover-icon');

      if (isFavoriteByIndex(cell.tableRowIndex!)) {
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
          ...showActions ? [unfavoriteIcon, addToFavorites, editIcon, deleteIcon]: [],
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
};

export const setGridColumnsRendering = (grid: DG.Grid) => {
  const actionsCol = grid.columns.byName(ACTIONS_COLUMN_NAME);
  if (actionsCol) {
    actionsCol.cellType = 'html';
    actionsCol.width = 45;
  }

  const favCol = grid.columns.byName(FAVORITE_COLUMN_NAME);
  if (favCol) {
    favCol.cellType = 'html';
    favCol.width = 30;
  }

  const expCol = grid.columns.byName(EXP_COLUMN_NAME);
  if (expCol) {
    expCol.cellType = 'html';
    expCol.width = 30;
  }

  const tagsColumn = grid.columns.byName(TAGS_COLUMN_NAME);
  if (tagsColumn) {
    tagsColumn.cellType = 'html';
    tagsColumn.width = 90;
  }

  const startedCol = grid.columns.byName(STARTED_COLUMN_NAME);
  if (startedCol)
    startedCol.width = 110;

  const idCol = grid.columns.byName(ID_COLUMN_NAME);
  if (idCol)
    idCol.cellType = 'html';
};

const camel2title = (camelCase: string) => camelCase
  .replace(/([A-Z])/g, (match) => ` ${match.toLowerCase()}`)
  .trim()
  .replace(/^./, (match) => match.toUpperCase());

export const getColumnName = (key: string) => {
  return camel2title(key);
};

export const extractStringValue = (run: DG.FuncCall, key: string) => {
  if (key === AUTHOR_COLUMN_NAME) return run.author?.friendlyName ?? ' ';

  const val =
    (run as any)[key] ??
    run.inputs[key] ??
    run.outputs[key] ??
    run.options[key] ??
    null;

  return val?.toString() ?? '';
};

export const getRunsDfFromList = async (
  runs: Map<string, DG.FuncCall>,
  func: DG.Func,
  options?: HistoryOptions,
) => {
  const newRuns = [...runs.values()];

  const getColumnByProp = (prop: DG.Property) => {
    if (prop.propertyType === DG.TYPE.DATE_TIME) {
      return DG.Column.dateTime(prop.caption ?? getColumnName(prop.name), newRuns.length)
        // Workaround for https://reddata.atlassian.net/browse/GROK-15286
        .init((idx) => (<any>window).grok_DayJs_To_DateTime(newRuns[idx].inputs[prop.name]));
    }

    return DG.Column.fromType(
      prop.propertyType as any,
      prop.caption ?? getColumnName(prop.name),
      newRuns.length,
    ).init((idx) =>
      (options?.propFuncs?.[prop.name])?.(newRuns[idx]) ??
      extractStringValue(newRuns[idx], prop.name),
    );
  };

  const getColumnByName = (key: string) => {
    if (key === STARTED_COLUMN_NAME) {
      return DG.Column.dateTime(getColumnName(key), newRuns.length)
        // Workaround for https://reddata.atlassian.net/browse/GROK-15286
        .init((idx) =>
          (<any>window).grok_DayJs_To_DateTime(getStartedOrNull(newRuns[idx]) ?
            newRuns[idx].started.utc(true): dayjs.unix(newRuns[idx].options['createdOn'])),
        );
    }

    if (key === COMPLETE_COLUMN_NAME) {
      return DG.Column.bool(getColumnName(key), newRuns.length)
        // Workaround for https://reddata.atlassian.net/browse/GROK-15286
        .init((idx) => !isIncomplete(newRuns[idx]));
    }

    return DG.Column.fromStrings(
      getColumnName(key),
      newRuns.map((run) => (options?.propFuncs?.[key])?.(run) ?? extractStringValue(run, key)),
    );
  };

  const favoritesRecord: Record<string, string> = await grok.dapi.userDataStorage.get(getFavStorageName(func)) ?? {};
  const favorites = Object.keys(favoritesRecord);

  const getColumn = (key: string) => {
    const prop =
    func.inputs.find((prop) => prop.name === key) ??
    func.outputs.find((prop) => prop.name === key);
    if (prop)
      return getColumnByProp(prop);
    else
      return getColumnByName(key);
  };

  const newRunsGridDf = DG.DataFrame.fromColumns([
    DG.Column.string(EXP_COLUMN_NAME, newRuns.length).init((idx) => {
      const immutableTags = newRuns[idx].options['immutable_tags'] as string[];
      return immutableTags && immutableTags.includes(EXPERIMENTAL_TAG) ? 'Experimental': 'Simulated';
    }),
    ...options?.isHistory ?
      [DG.Column.bool(FAVORITE_COLUMN_NAME, newRuns.length)
        .init((idx) => favorites.includes(newRuns[idx].id))]: [],
    ...options?.showActions ? [DG.Column.string(ACTIONS_COLUMN_NAME, newRuns.length).init('')]: [],
    getColumnByName(STARTED_COLUMN_NAME),
    getColumnByName(COMPLETE_COLUMN_NAME),
    getColumnByName(AUTHOR_COLUMN_NAME),
    DG.Column.string(TAGS_COLUMN_NAME, newRuns.length).init((idx) =>
      newRuns[idx].options['tags'] ? newRuns[idx].options['tags'].join(','): '',
    ),
    DG.Column.string(TITLE_COLUMN_NAME, newRuns.length).init((idx) => newRuns[idx].options['title']),
    DG.Column.string(DESC_COLUMN_NAME, newRuns.length).init((idx) => newRuns[idx].options['description']),
    DG.Column.string(VERSION_COLUMN_NAME, newRuns.length).init((idx) => newRuns[idx].options['version']),
  ]);

  getVisibleProps(func).map((key) => getColumn(key)).forEach((col) => {
    col.name = newRunsGridDf.columns.getUnusedName(col.name);
    newRunsGridDf.columns.add(col, false);
  });

  newRunsGridDf.columns.add(DG.Column.fromStrings(ID_COLUMN_NAME, newRuns.map((newRun) => newRun.id)));

  if (!options?.allowOtherVersions && options?.version) {
    const rowMask = DG.BitSet.create(newRunsGridDf.rowCount, () => false);
    for (let idx = 0; idx < newRunsGridDf.rowCount; idx++) {
      if (options?.version == newRunsGridDf.get(VERSION_COLUMN_NAME, idx))
        rowMask.set(idx, true);
    }
    return newRunsGridDf.clone(rowMask);
  }
  return newRunsGridDf;
};
