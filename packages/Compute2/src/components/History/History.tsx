import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {IconFA, ifOverlapping, ToggleInput, Viewer} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {ID_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-components/src/history-input';
import {FAVORITE_COLUMN_NAME, COMPLETE_COLUMN_NAME, STARTED_COLUMN_NAME, AUTHOR_COLUMN_NAME, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME, DESC_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {HistoricalRunEdit, HistoricalRunsDelete} from '@datagrok-libraries/compute-utils/shared-components/src/history-dialogs';
import {filter, mergeMap, take, toArray} from 'rxjs/operators';
import {from} from 'rxjs';
import wu from 'wu';
import {watchExtractedObservable} from '@vueuse/rxjs';

const GRID_INITED_EVENT = 'd4-grid-initialized';

export const History = Vue.defineComponent({
  name: 'History',
  props: {
    func: {
      type: DG.Func,
      required: true,
    },
    fallbackText: {
      type: String,
      default: 'No historical runs found',
    },
    showActions: {
      type: Boolean,
      default: false,
    },
    showBatchActions: {
      type: Boolean,
      default: false,
    },
    isHistory: {
      type: Boolean,
      default: false,
    },
  },
  emits: {
    runChosen: (_chosenCall: DG.FuncCall) => true,
    compare: (_ids: string[]) => true,
    afterRunEdited: (_editedCall: DG.FuncCall) => true,
    afterRunDeleted: (_deletedCall: DG.FuncCall) => true,
  },
  setup(props, {emit}) {
    const isLoading = Vue.ref(true);

    const showFilters = Vue.ref(false);
    const showInputs = Vue.ref(false);
    const showMetadata = Vue.ref(true);
    const hasSelected = Vue.ref(false);

    const historicalRuns = Vue.shallowRef(new Map<string, DG.FuncCall>);
    const func = Vue.computed(() => Vue.markRaw(props.func));

    const currentSelection = new Set<DG.FuncCall>();

    const refresh = async () => {
      isLoading.value = true;
      try {
        const newHistoricalRuns = await historyUtils.pullRunsByName(
          func.value.name,
          [{author: grok.shell.user}],
          {},
          ['session.user', 'options'],
        );
        historicalRuns.value.clear();

        newHistoricalRuns.reduce((acc, run) => {
          acc.set(run.id, run);
          return acc;
        }, historicalRuns.value);
        Vue.triggerRef(historicalRuns);
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isLoading.value = false;
      }
    };

    Vue.watch(func, () => refresh(), {immediate: true});

    const defaultDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(ID_COLUMN_NAME, []),
      DG.Column.bool(COMPLETE_COLUMN_NAME, 0),
      DG.Column.dateTime(STARTED_COLUMN_NAME, 0),
      DG.Column.string(AUTHOR_COLUMN_NAME, 0),
      DG.Column.string(TAGS_COLUMN_NAME, 0),
      DG.Column.string(TITLE_COLUMN_NAME, 0),
      DG.Column.string(DESC_COLUMN_NAME, 0),
    ]);

    const getRunByIdx = (idx: number) => {
      if (idx < 0) return;

      const run = historicalRuns.value.get(historicalRunsDf.value.get(ID_COLUMN_NAME, idx));
      return run ? Vue.markRaw(run) : undefined;
    };

    const isFavoriteByIndex = (idx: number) => {
      return historicalRunsDf.value.get(FAVORITE_COLUMN_NAME, idx);
    };

    const updateRun = (updatedRun: DG.FuncCall) => {
      Vue.markRaw(updatedRun);
      historicalRuns.value.set(updatedRun.id, updatedRun);
      Vue.triggerRef(historicalRuns);
    };

    const showEditDialog = async (funcCall: DG.FuncCall, isFavorite: boolean) => {
      const editDialog = new HistoricalRunEdit(funcCall, isFavorite);
      editDialog.show({center: true, width: 500});
      const editOptions = await editDialog.onMetadataEdit.pipe(take(1)).toPromise();
      try {
        if (!props.isHistory)
          return updateRun(funcCall);
        if (editOptions.favorite !== 'same')
          await Utils.saveIsFavorite(funcCall, (editOptions.favorite === 'favorited'));
        const fullCall = await historyUtils.loadRun(funcCall.id, false, false);
        if (editOptions.title) fullCall.options['title'] = editOptions.title;
        if (editOptions.description) fullCall.options['description'] = editOptions.description;
        if (editOptions.tags) fullCall.options['tags'] = editOptions.tags;
        await historyUtils.saveRun(fullCall);
        updateRun(fullCall);
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const onEditClick = (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      showEditDialog(
        run,
        isFavoriteByIndex(cell.tableRowIndex!),
      );
    };

    const onFavoriteClick = async (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      await Utils.saveIsFavorite(run, true);
      updateRun(run);
    };

    const onUnfavoriteClick = async (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      await Utils.saveIsFavorite(run, false);
      updateRun(run);
    };

    const deleteRun = async (id: string) => {
      try {
        const loadedRun = await historyUtils.loadRun(id, true, false);
        if (props.isHistory)
          await historyUtils.deleteRun(loadedRun);
        historicalRuns.value.delete(id);
        Vue.triggerRef(historicalRuns);
        await Utils.saveIsFavorite(loadedRun, false);
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const onDeleteClick = async (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      const setToDelete = new Set([run]);
      const deleteDialog = new HistoricalRunsDelete(setToDelete);
      deleteDialog.show({center: true, width: 500});

      await deleteDialog.onFuncCallDelete.pipe(take(1)).toPromise();
      try {
        await Promise.all(
          wu(setToDelete.values()).map(async (funcCall) => {
            await deleteRun(funcCall.id);
          }));
      } catch (e: any) {
        grok.shell.error(e);
      }
    };

    const historicalRunsDf = Vue.shallowRef(Vue.markRaw(defaultDf));

    Vue.watch([historicalRuns, func], async ([historicalRuns, func]) => {
      const df = await Utils.getRunsDfFromList(
        historicalRuns,
        func,
        Vue.toValue(() => props),
      );
      historicalRunsDf.value = Vue.markRaw(df);
    });

    watchExtractedObservable(historicalRunsDf, (p) => p.onCurrentRowChanged, async () => {
      historicalRunsDf.value.rows.select(() => false);
      const chosenRun = getRunByIdx(historicalRunsDf.value.currentRowIdx);
      if (chosenRun) {
        const run = await historyUtils.loadRun(chosenRun.id, false, false);
        emit('runChosen', run ? Vue.markRaw(run) : run);
      };
    });

    watchExtractedObservable(historicalRunsDf, (p) => p.onSelectionChanged, () => {
      const selectedCalls: DG.FuncCall[] = [...historicalRunsDf.value.selection
        .getSelectedIndexes()]
        .map((idx) => getRunByIdx(idx)!);
      currentSelection.clear();
      for (const call of selectedCalls)
        currentSelection.add(call);
      hasSelected.value = selectedCalls.length > 0;
    });

    const currentGrid = Vue.shallowRef<null | DG.Grid>(null);

    const updateVisibleColumns = () => {
      if (!currentGrid.value) return;

      const tagCol = currentGrid.value.dataFrame.getCol(TAGS_COLUMN_NAME);
      const cols = [
        ID_COLUMN_NAME,
        ...showMetadata.value ? [STARTED_COLUMN_NAME, COMPLETE_COLUMN_NAME]: [],
        ...showMetadata.value && tagCol.stats.missingValueCount < tagCol.length ? [TAGS_COLUMN_NAME]: [],
        ...showMetadata.value ? [TITLE_COLUMN_NAME]: [],
        ...showMetadata.value ? [DESC_COLUMN_NAME]: [],
        ...showInputs.value && currentFunc.value ? Utils.getVisibleProps(currentFunc.value)
          .map((key) => {
            const param = currentFunc.value.inputs.find((prop) => prop.name === key) ??
              currentFunc.value.outputs.find((prop) => prop.name === key);

            if (param)
              return param.caption ?? Utils.getColumnName(param.name);
            else
              return Utils.getColumnName(key);
          }): [],
      ];
      currentGrid.value.columns.setVisible(cols);
      currentGrid.value.columns.setOrder(cols);
    };

    Vue.watch([showMetadata, showInputs, historicalRunsDf, currentGrid], () => updateVisibleColumns(), {flush: 'post'});

    const handleGridRendering = async (grid?: DG.Grid) => {
      if (!grid) return;

      await grok.events.onEvent(GRID_INITED_EVENT).pipe(
        filter((initedGrid) => initedGrid === grid),
        take(1),
      ).toPromise();

      currentGrid.value = Vue.markRaw(grid);
      grid.sort([STARTED_COLUMN_NAME], [false]);
      grid.props.rowHeight = 46;

      for (let i = 0; i < grid.columns.length; i++) {
        const col = grid.columns.byIndex(i);
        if (col && col.column?.type === DG.TYPE.DATE_TIME)
          col.format = 'MMM d, h:mm tt';
      }

      Utils.setGridColumnsRendering(grid);

      Utils.setGridCellRendering(
        grid,
        historicalRuns.value,
        onEditClick,
        onDeleteClick,
        onFavoriteClick,
        onUnfavoriteClick,
        true,
      );
    };

    const fallbackText = (
      <div class='p-1'>
        {props.fallbackText}
        <IconFA
          name='sync'
          tooltip={'Refresh'}
          onClick={refresh}
          style={{alignContent: 'center', paddingLeft: '10px'}}
        />
      </div>
    );

    const deleteSelected = async () => {
      if (isLoading.value || currentSelection.size === 0)
        return;
      const deleteDialog = new HistoricalRunsDelete(currentSelection);
      deleteDialog.show({center: true, width: 500});
      await deleteDialog.onFuncCallDelete.pipe(take(1)).toPromise();
      try {
        isLoading.value = true;
        await from(currentSelection).pipe(
          mergeMap(call => from(deleteRun(call.id)), 5),
          toArray(),
        ).toPromise();
      } finally {
        isLoading.value = false;
      }
    };

    const currentFunc = Vue.computed(() => Vue.markRaw(props.func));
    const isHistory = Vue.computed(() => props.isHistory);
    const visibleFilterColumns = Vue.computed(() => {
      const currentDf = historicalRunsDf.value;
      const tagCol = currentDf.getCol(TAGS_COLUMN_NAME);

      return [
        ...showMetadata.value &&
        (currentDf.col(FAVORITE_COLUMN_NAME)?.categories.length ?? 0) > 1 ?
          [FAVORITE_COLUMN_NAME]: [],
        ...showMetadata.value ? [STARTED_COLUMN_NAME, COMPLETE_COLUMN_NAME]:[],
        ...showMetadata.value &&
        currentDf.getCol(AUTHOR_COLUMN_NAME).categories.length > 1 ? [AUTHOR_COLUMN_NAME]: [],
        ...showMetadata.value &&
        tagCol.stats.missingValueCount < tagCol.length ? [TAGS_COLUMN_NAME]: [],
        ...showMetadata.value &&
        currentDf.getCol(TITLE_COLUMN_NAME).categories.length > 1 ? [TITLE_COLUMN_NAME]: [],
        ...showMetadata.value &&
        currentDf.getCol(DESC_COLUMN_NAME).categories.length > 1 ? [DESC_COLUMN_NAME]: [],
        ...currentFunc.value && showInputs.value ? Utils.getVisibleProps(currentFunc.value)
          .map((key) => {
            const param = currentFunc.value.inputs.find((prop) => prop.name === key) ??
            currentFunc.value.outputs.find((prop) => prop.name === key);

            if (param)
              return param.caption ?? Utils.getColumnName(param.name);
            else
              return Utils.getColumnName(key);
          })
          .filter((columnName) => {
            return !isHistory.value ||
            currentDf.getCol(columnName).categories.length > 1;
          })
          .map((columnName) => columnName): [],
      ];
    });

    return () => {
      const controls = <div style={{display: 'flex', justifyContent: 'space-between', padding: '0px 6px'}}>
        <div style={{'display': 'flex', 'padding': '6px 0px', 'gap': '6px'}}>
          <IconFA
            name='filter'
            tooltip={showFilters.value ? 'Hide filters': 'Show filters'}
            faStyle={showFilters.value ? 'fal': 'fad'}
            onClick={() => showFilters.value = !showFilters.value}
            style={{alignContent: 'center'}}
          />
          <ToggleInput
            caption='Metadata'
            value={showMetadata.value}
            onUpdate:value={(val) => showMetadata.value = val}
          />
          <ToggleInput
            caption='Params'
            value={showInputs.value}
            onUpdate:value={(val) => showInputs.value = val}
          />
        </div>
        <div style={{display: 'flex', 'padding': '6px 0px', 'gap': '6px'}}>
          <IconFA
            name='sync'
            tooltip={'Refresh'}
            onClick={refresh}
            style={{alignContent: 'center'}}
          />
          { hasSelected.value && <IconFA
            name='trash-alt'
            tooltip={'Delete'}
            onClick={deleteSelected}
            style={{alignContent: 'center'}}
          />}
        </div>
      </div>;
      const grid = <Viewer
        type='Grid'
        dataFrame={historicalRunsDf.value}
        style={{height: '100%', width: '100%', minHeight: '300px'}}
        onViewerChanged={(viewer) => handleGridRendering(viewer as DG.Grid | undefined)}
        options={{
          'showCurrentRowIndicator': true,
          'showCurrentCellOutline': false,
          'allowEdit': false,
          'allowBlockSelection': false,
          'showRowHeader': false,
          'showColumnLabels': true,
          'extendLastColumn': false,
        }}
      />;
      const filters = <Viewer
        type='Filters'
        dataFrame={historicalRunsDf.value}
        style={{
          flex: '1', width: '100%', overflow: 'hidden',
        }}
        options={{
          'columnNames': visibleFilterColumns.value,
          'showHeader': false,
          'allowEdit': false,
        }}
      />;

      return Vue.withDirectives(<div style={{overflow: 'hidden', height: '100%'}}>
        { historicalRuns.value.size === 0 ?
          fallbackText:
          <div style={{
            display: 'flex',
            flexDirection: 'column',
            width: '100%',
            height: '100%',
          }}>
            <div style={{display: 'flex', flexDirection: 'column', flex: '1'}}>
              { controls }
              { grid }
            </div>
            { showFilters.value && filters }
          </div>
        }
      </div>, [[ifOverlapping, isLoading.value, 'Loading previous runs...']]);
    };
  },
});
