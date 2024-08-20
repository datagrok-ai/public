import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computed, defineComponent, PropType, ref, shallowReactive, shallowRef, ShallowRef, toValue, watch, watchEffect} from 'vue';
import {IconFA, ToggleInput, Viewer} from '@datagrok-libraries/webcomponents-vue/src';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getRunsDfFromList, styleHistoryFilters, styleHistoryGrid} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {computedAsync} from '@vueuse/core';
import {ID_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-components/src/history-input';
import {EXP_COLUMN_NAME, FAVORITE_COLUMN_NAME, ACTIONS_COLUMN_NAME, COMPLETE_COLUMN_NAME, STARTED_COLUMN_NAME, AUTHOR_COLUMN_NAME, TAGS_COLUMN_NAME, TITLE_COLUMN_NAME, DESC_COLUMN_NAME} from '@datagrok-libraries/compute-utils/shared-utils/consts';

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
    onRunChosen: (id: string) => id,
    onComparison: (ids: string[]) => ids,
    afterRunEdited: (editedCall: DG.FuncCall) => editedCall,
    afterRunDeleted: (deletedCall: DG.FuncCall) => deletedCall,
  },
  setup(props) {
    const isLoading = ref(true);

    const isFullMode = ref(false);
    const showFilters = ref(true);
    const showInputs = ref(true);
    const showMetadata = ref(true);

    const historicalRuns = shallowReactive(new Map<string, DG.FuncCall>);

    watch(() => props.func.name, () => {
      isLoading.value = true;

      historyUtils.pullRunsByName(props.func.name, [{author: grok.shell.user}], {}, ['session.user', 'options'])
        .then((newHistoricalRuns) => {
          historicalRuns.clear();

          newHistoricalRuns.reduce((acc, run) => {
            acc.set(run.id, run);
            return acc;
          }, historicalRuns);
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

    watchEffect(() => {
      if (!currentGrid.value || !currentFilters.value) return;

      const func = historicalRuns.values().next().value.func as DG.Func | undefined;

      styleHistoryGrid(
        currentGrid.value, 
        isFullMode.value,
        showInputs.value,
        showMetadata.value,
        func,
      );

      styleHistoryFilters(
        currentFilters.value,
        showMetadata.value,
        showInputs.value,
        props.isHistory,
        func,
      );
    });
    
    const historicalRunsDf = computedAsync(async () => {
      return await getRunsDfFromList(
        historicalRuns, 
        toValue(() => props),
      );
    }, defaultDf);

    return () => {
      const controls = <div style={{display: 'flex', justifyContent: 'space-between', padding: '0px 6px'}}>
        <div style={{'display': 'flex', 'padding': '6px 0px', 'gap': '6px'}}>
          <IconFA 
            name={isFullMode.value ? 'compress-alt': 'expand-alt'} 
            tooltip={isFullMode.value ? 'Switch to compact mode': 'Switch to full mode'}
            onClick={() => isFullMode.value = !isFullMode.value}
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
        style={{height: '100%', width: '100%'}}
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
            flexDirection: isFullMode.value ? 'row': 'column', 
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
