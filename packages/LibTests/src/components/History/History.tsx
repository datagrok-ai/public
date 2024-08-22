import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, nextTick, onMounted, PropType, ref, shallowReactive, shallowRef, ShallowRef, toValue, triggerRef, watch, watchEffect} from 'vue';
import {IconFA, ToggleInput, Viewer} from '@datagrok-libraries/webcomponents-vue/src';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import * as Utils from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {computedAsync} from '@vueuse/core';
import {ID_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-components/src/history-input';
import {EXP_COLUMN_NAME, FAVORITE_COLUMN_NAME, ACTIONS_COLUMN_NAME, COMPLETE_COLUMN_NAME, STARTED_COLUMN_NAME, AUTHOR_COLUMN_NAME, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME, DESC_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-utils/consts';
import {HistoricalRunEdit, HistoricalRunsDelete} from '@datagrok-libraries/compute-utils/shared-components/src/history-dialogs';
import {take} from 'rxjs/operators';
import wu from 'wu';
import {useObservable, useSubscription, watchExtractedObservable} from '@vueuse/rxjs';

export const History = defineComponent({
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
    propFuncs: {
      type: Object as PropType<Record<string, (currentRun: DG.FuncCall) => string>>,
      default: {},
    }, 
  },
  emits: {
    runChosen: (chosenCall: DG.FuncCall) => chosenCall,
    compare: (ids: string[]) => ids,
    afterRunEdited: (editedCall: DG.FuncCall) => editedCall,
    afterRunDeleted: (deletedCall: DG.FuncCall) => deletedCall,
  },
  setup(props, {emit}) {
    const isLoading = ref(true);

    const isCompactMode = ref(true);
    const showFilters = ref(true);
    const showInputs = ref(true);
    const showMetadata = ref(true);

    const historicalRuns = shallowRef(new Map<string, DG.FuncCall>);

    watch(() => props.func.name, () => {
      isLoading.value = true;

      historyUtils.pullRunsByName(props.func.name, [{author: grok.shell.user}], {}, ['session.user', 'options'])
        .then((newHistoricalRuns) => {
          historicalRuns.value.clear();

          newHistoricalRuns.reduce((acc, run) => {
            acc.set(run.id, run);
            return acc;
          }, historicalRuns.value);
          triggerRef(historicalRuns);
        })
        .catch((e) => grok.shell.error(e))
        .finally(() => isLoading.value = false);
    }, {immediate: true});

    const defaultDf = DG.DataFrame.fromColumns([
      DG.Column.bool(EXP_COLUMN_NAME, 0),
      ...props.isHistory ? [DG.Column.bool(FAVORITE_COLUMN_NAME, 0)]: [],
      ...props.showActions ? [DG.Column.string(ACTIONS_COLUMN_NAME, 0)]: [],
      DG.Column.bool(COMPLETE_COLUMN_NAME, 0),
      DG.Column.dateTime(STARTED_COLUMN_NAME, 0),
      DG.Column.string(AUTHOR_COLUMN_NAME, 0),
      DG.Column.string(TAGS_COLUMN_NAME, 0),
      DG.Column.string(TITLE_COLUMN_NAME, 0),
      DG.Column.string(DESC_COLUMN_NAME, 0),
      DG.Column.fromStrings(ID_COLUMN_NAME, []),
    ]);

    const currentGrid = shallowRef(null as null | DG.Grid);
    const currentFilters = shallowRef(null as null | DG.FilterGroup);

    const getRunByIdx = (idx: number) => {
      if (idx < 0) return;

      return historicalRuns.value.get(historicalRunsDf.value.get(ID_COLUMN_NAME, idx));
    };

    const isFavoriteByIndex = (idx: number) => {
      return historicalRunsDf.value.get(FAVORITE_COLUMN_NAME, idx);
    };

    const updateRun = (updatedRun: DG.FuncCall) => {
      historicalRuns.value.set(updatedRun.id, updatedRun);
      triggerRef(historicalRuns);
    };

    const showEditDialog = (funcCall: DG.FuncCall, isFavorite: boolean) => {
      const editDialog = new HistoricalRunEdit(funcCall, isFavorite);
  
      editDialog.onMetadataEdit.pipe(take(1)).subscribe(async (editOptions) => {
        if (!props.isHistory)
          updateRun(funcCall);
        else {
          return ((editOptions.favorite !== 'same') ?
            Utils.saveIsFavorite(funcCall, (editOptions.favorite === 'favorited')) :
            Promise.resolve())
            .then(() => historyUtils.loadRun(funcCall.id, false))
            .then((fullCall) => {
              if (editOptions.title) fullCall.options['title'] = editOptions.title;
              if (editOptions.description) fullCall.options['description'] = editOptions.description;
              if (editOptions.tags) fullCall.options['tags'] = editOptions.tags;
  
              return [historyUtils.saveRun(fullCall), fullCall] as const;
            })
            .then(([, fullCall]) => {
              updateRun(fullCall);
            })
            .catch((err) => {
              grok.shell.error(err);
            });
        }
      });
      editDialog.show({center: true, width: 500});
    };
    
    const onEditClick = (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      showEditDialog(
        run,
        isFavoriteByIndex(cell.tableRowIndex!),
      );
    };

    const onFavoriteClick = (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      Utils.saveIsFavorite(run, true).then(() => updateRun(run));
    };

    const onUnfavoriteClick = (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      Utils.saveIsFavorite(run, false).then(() => updateRun(run));
    };

    const deleteRun = async (id: string) => {
      return historyUtils.loadRun(id, true)
        .then(async (loadedRun) => {
          return [
            await (props.isHistory ? historyUtils.deleteRun(loadedRun): Promise.resolve()),
            loadedRun,
          ] as const;
        })
        .then(([, loadedRun]) => {
          historicalRuns.value.delete(id);
          triggerRef(historicalRuns);

          return loadedRun;
        })
        .then((loadedRun) => {
          Utils.saveIsFavorite(loadedRun, false);
        })
        .catch((e) => {
          grok.shell.error(e);
        });
    };

    const onDeleteClick = (cell: DG.GridCell) => {
      const run = getRunByIdx(cell.tableRowIndex!)!;
      const setToDelete = new Set([run]);
      const deleteDialog = new HistoricalRunsDelete(setToDelete);

      deleteDialog.onFuncCallDelete.pipe(
        take(1),
      ).subscribe(async () => {
        // ui.setUpdateIndicator(this.root, true);
        try {
          await Promise.all(
            wu(setToDelete.values()).map(async (funcCall) => {
              await deleteRun(funcCall.id);

              return Promise.resolve();
            }));
        } catch (e: any) {
          grok.shell.error(e);
        } finally {
          // ui.setUpdateIndicator(this.root, false);
        }
      });
      deleteDialog.show({center: true, width: 500});
    };

    watch(currentGrid, () => {
      Utils.setGridCellRendering(
        currentGrid.value!,
        historicalRuns.value,
        onEditClick,
        onDeleteClick,
        onFavoriteClick,
        onUnfavoriteClick,
        true,
      );
    });

    const historicalRunsDf = computedAsync(async () => {
      const df = await Utils.getRunsDfFromList(
        historicalRuns.value, 
        props.func,
        toValue(() => props),
      );
      return df;
    }, defaultDf);

    watchExtractedObservable(historicalRunsDf, (p) => p.onCurrentRowChanged, async () => {
      historicalRunsDf.value.rows.select(() => false);
      const chosenRun = getRunByIdx(historicalRunsDf.value.currentRowIdx);
      if (chosenRun) emit('runChosen', await historyUtils.loadRun(chosenRun.id));
    });
  
    const applyStyles = () => {
      const func = historicalRuns.value.values().next().value?.func as DG.Func | undefined;

      Utils.styleHistoryGrid(
        currentGrid.value!, 
        isCompactMode.value,
        showInputs.value,
        showMetadata.value,
        func,
      );

      Utils.styleHistoryFilters(
        currentFilters.value!,
        showMetadata.value,
        showInputs.value,
        props.isHistory,
        func,
      );
    };

    watch([showInputs, showMetadata, isCompactMode], () => {
      applyStyles();
    });

    watch([historicalRuns, currentGrid, currentFilters], () => {
      if (!currentGrid.value || !currentFilters.value) return;

      setTimeout(() => {
        applyStyles();
      }, 100);
    });

    return () => {
      const controls = <div style={{display: 'flex', justifyContent: 'space-between', padding: '0px 6px'}}>
        <div style={{'display': 'flex', 'padding': '6px 0px', 'gap': '6px'}}>
          <IconFA 
            name={isCompactMode.value ? 'expand-alt': 'compress-alt'} 
            tooltip={isCompactMode.value ? 'Switch to full mode': 'Switch to compact mode'}
            onClick={() => isCompactMode.value = !isCompactMode.value}
            style={{alignContent: 'center'}}
          />
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
        <div style={{display: 'flex'}}>
          
        </div>
      </div>;
      const grid = <Viewer 
        type='Grid'
        dataFrame={historicalRunsDf.value} 
        style={{height: '100%', width: '100%', minHeight: '300px'}}
        onViewerChanged={(viewer) => currentGrid.value = viewer as DG.Grid}
      />;
      const filters = <Viewer 
        type='Filters' 
        dataFrame={historicalRunsDf.value} 
        style={{height: '100%', width: '100%', display: !showFilters.value ? 'none': 'block'}}
        onViewerChanged={(viewer) => currentFilters.value = viewer as DG.FilterGroup}
      />;

      return <div style={{height: '100%'}}>
        {isLoading.value ? 
          <span> Loading... </span>: 
          <div style={{
            display: 'flex', 
            flexDirection: isCompactMode.value ? 'column': 'row', 
            height: '100%', width: '100%',
          }}> 
            <div style={{display: 'flex', flexDirection: 'column', height: '100%', width: '100%'}}>
              { controls }
              { grid }
            </div>
            { filters }
          </div>
        }
      </div>;
    };
  },
});
