import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {IconFA, Viewer, ResizeHandle, ToggleInput, ifOverlapping} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';

import {History} from '../History/History';
import {
  ComparisonTarget, MatchConfidence, matchScalarTargets, matchColumnTargets, getEntryStatuses,
  normalizeName, nameSimilarity, FUZZY_NAME_THRESHOLD,
} from './comparison-core';
import {
  ComparisonEntry, ComparisonMode, EntrySourceKind, RUN_COLUMN,
  entryFromFuncCall, entryFromDataFrame, buildScalarComparison, buildColumnComparison,
} from './data-extraction';
import './RunComparison.css';

const SOURCE_BADGES: Record<EntrySourceKind, {label: string, color: string}> = {
  workflow: {label: 'workflow', color: '#7b6fb3'},
  function: {label: 'function', color: '#4a90d9'},
  raw: {label: 'raw', color: '#8a8a8a'},
};

const CONFIDENCE_COLORS: Record<MatchConfidence, string> = {
  exact: '#3cb173',
  normalized: '#d9a544',
  fuzzy: '#d97b44',
};

export const RunComparison = Vue.defineComponent({
  name: 'RunComparison',
  props: {
    roleOnlyFilter: {
      type: Boolean,
      default: false,
    },
  },
  setup(props) {
    const isLoadingModels = Vue.ref(true);
    const isAddingRuns = Vue.ref(false);
    const models = Vue.shallowRef<DG.Func[]>([]);
    const selectedModel = Vue.shallowRef<DG.Func | undefined>(undefined);
    const modelSearch = Vue.ref('');
    const modelDropdownOpen = Vue.ref(false);
    const modelDropdownRef = Vue.ref<HTMLElement>();

    const allowFloatIndex = Vue.ref(false);
    const allowDatetimeIndex = Vue.ref(false);

    const historySelection = Vue.shallowRef<DG.FuncCall[]>([]);
    const entries = Vue.shallowRef<ComparisonEntry[]>([]);
    const baselineId = Vue.ref('');
    // entryId -> tablePath -> index column name
    const indexSelection = Vue.ref<Record<string, Record<string, string>>>({});
    const selectedTargetKey = Vue.ref('');
    const mode = Vue.ref<ComparisonMode>('values');

    const historyHeight = Vue.ref(320);
    const sidebarWidth = Vue.ref(360);

    Vue.onMounted(async () => {
      try {
        const list = await grok.dapi.functions.filter(getModelFilter(props.roleOnlyFilter)).list();
        models.value = list.map((f) => Vue.markRaw(f));
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isLoadingModels.value = false;
      }
    });

    const onClickOutside = (e: MouseEvent) => {
      if (modelDropdownRef.value && !modelDropdownRef.value.contains(e.target as Node))
        modelDropdownOpen.value = false;
    };
    Vue.onMounted(() => document.addEventListener('mousedown', onClickOutside));
    Vue.onUnmounted(() => document.removeEventListener('mousedown', onClickOutside));

    const modelLabel = (f: DG.Func) => f.friendlyName ?? f.name;
    const filteredModels = Vue.computed(() => {
      const query = modelSearch.value.toLowerCase();
      return models.value.filter((f) => !query || modelLabel(f).toLowerCase().includes(query));
    });

    const tableInput = ui.input.table('Add table');
    tableInput.onChanged.subscribe(() => {
      const df = tableInput.value;
      if (!df)
        return;
      addEntry(entryFromDataFrame(df));
      tableInput.value = null;
    });

    const addEntry = (entry: ComparisonEntry) => {
      if (entries.value.some((existing) => existing.id === entry.id))
        return;
      entries.value = [...entries.value, Vue.markRaw(entry)];
      if (!baselineId.value)
        baselineId.value = entry.id;
    };

    const removeEntry = (id: string) => {
      entries.value = entries.value.filter((entry) => entry.id !== id);
      if (baselineId.value === id)
        baselineId.value = entries.value[0]?.id ?? '';
      const {[id]: _removed, ...rest} = indexSelection.value;
      indexSelection.value = rest;
    };

    const addSelectedRuns = async () => {
      if (isAddingRuns.value || historySelection.value.length === 0)
        return;
      isAddingRuns.value = true;
      try {
        for (const call of historySelection.value) {
          if (entries.value.some((entry) => entry.id === call.id))
            continue;
          const fullCall = await historyUtils.loadRun(call.id);
          addEntry(await entryFromFuncCall(fullCall));
        }
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isAddingRuns.value = false;
      }
    };

    const setIndexColumn = (entryId: string, tablePath: string, columnName: string) => {
      indexSelection.value = {
        ...indexSelection.value,
        [entryId]: {...indexSelection.value[entryId], [tablePath]: columnName},
      };
    };

    // float/datetime indexes are edge cases prone to alignment noise, so they are opt-in
    const isAllowedIndexType = (type?: string) => {
      if (type === DG.COLUMN_TYPE.INT || type === DG.COLUMN_TYPE.BIG_INT || type === DG.COLUMN_TYPE.STRING)
        return true;
      if (type === DG.COLUMN_TYPE.DATE_TIME)
        return allowDatetimeIndex.value;
      if (type === DG.COLUMN_TYPE.FLOAT)
        return allowFloatIndex.value;
      return false;
    };

    const indexColumnType = (entryId: string, tablePath: string, columnName: string) => {
      const entry = entries.value.find((item) => item.id === entryId);
      const table = entry?.nodes.tables.find((item) => item.path === tablePath);
      return table?.columns.find((col) => col.name === columnName)?.type;
    };

    const indexColumnsMap = Vue.computed(() => {
      const map = new Map<string, Map<string, string>>();
      for (const [entryId, tables] of Object.entries(indexSelection.value)) {
        const tableMap = new Map<string, string>();
        for (const [tablePath, columnName] of Object.entries(tables)) {
          if (columnName && isAllowedIndexType(indexColumnType(entryId, tablePath, columnName)))
            tableMap.set(tablePath, columnName);
        }
        map.set(entryId, tableMap);
      }
      return map;
    });

    const targets = Vue.computed<ComparisonTarget[]>(() => {
      const nodes = entries.value.map((entry) => entry.nodes);
      return [
        ...matchScalarTargets(nodes),
        ...matchColumnTargets(nodes, indexColumnsMap.value),
      ].sort((a, b) => b.coverage - a.coverage);
    });

    const indexFilter = Vue.ref('');
    const targetFilter = Vue.ref('');

    // substring match on the displayed text, with fuzzy per-token fallback for typos
    const matchesFilter = (query: string, text: string) => {
      const q = normalizeName(query);
      if (!q)
        return true;
      const t = normalizeName(text);
      return t.includes(q) || t.split(' ').some((token) => nameSimilarity(token, q) >= FUZZY_NAME_THRESHOLD);
    };

    const filteredTargets = Vue.computed(() =>
      targets.value.filter((target) => matchesFilter(targetFilter.value, target.displayName)));

    const selectedTarget = Vue.computed(() =>
      targets.value.find((target) => target.key === selectedTargetKey.value) ?? null);

    const entryStatuses = Vue.computed(() => getEntryStatuses(
      entries.value.map((entry) => entry.nodes),
      selectedTarget.value,
      indexColumnsMap.value,
    ));

    const comparison = Vue.computed(() => {
      const target = selectedTarget.value;
      if (!target || entries.value.length < 2)
        return null;
      if (target.kind === 'scalar') {
        const result = buildScalarComparison(target, entries.value, baselineId.value);
        return Vue.markRaw({kind: 'scalar' as const, target, ...result});
      }
      const result = buildColumnComparison(target, entries.value, baselineId.value, mode.value);
      return result ? Vue.markRaw({kind: 'column' as const, target, ...result}) : null;
    });

    const scalarValueColumn = Vue.computed(() => {
      const target = selectedTarget.value;
      if (!target)
        return '';
      if (mode.value === 'delta')
        return 'Δ';
      if (mode.value === 'deltaPct')
        return 'Δ%';
      return target.displayName;
    });

    // tables that could participate in column comparison and need an index choice
    const indexTables = Vue.computed(() => entries.value.flatMap((entry) =>
      entry.nodes.tables.map((table) => {
        const candidates = table.columns.filter((col) => isAllowedIndexType(col.type));
        const stored = indexSelection.value[entry.id]?.[table.path] ?? '';
        return {
          entryId: entry.id,
          entryName: entry.name,
          table,
          candidates,
          current: candidates.some((col) => col.name === stored) ? stored : '',
        };
      }),
    ));

    const filteredIndexTables = Vue.computed(() => indexTables.value.filter((item) =>
      matchesFilter(indexFilter.value, `${item.entryName} · ${item.table.friendlyPath ?? item.table.path}`)));

    const suggestedIndex = (columns: {name: string, type: string}[]) =>
      columns.find((col) => col.type === DG.COLUMN_TYPE.DATE_TIME)?.name ??
      columns.find((col) => col.type === DG.COLUMN_TYPE.FLOAT || col.type === DG.COLUMN_TYPE.INT ||
        col.type === DG.COLUMN_TYPE.BIG_INT)?.name;

    const renderListFilter = (value: Vue.Ref<string>, placeholder: string) => (
      <input
        type='text'
        placeholder={placeholder}
        value={value.value}
        onInput={(e: Event) => value.value = (e.target as HTMLInputElement).value}
        style={{
          border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 6px',
          outline: 'none', width: '160px',
        }}
      />
    );

    const renderModelSelector = () => (
      <div ref={modelDropdownRef} style={{position: 'relative'}}>
        <div
          style={{
            padding: '4px 8px', cursor: 'pointer', border: '1px solid var(--grey-2)',
            borderRadius: '3px', userSelect: 'none', display: 'flex', justifyContent: 'space-between',
          }}
          onClick={() => { modelDropdownOpen.value = !modelDropdownOpen.value; modelSearch.value = ''; }}
        >
          <span>{selectedModel.value ? modelLabel(selectedModel.value) : 'Select model...'}</span>
          <span>&#9662;</span>
        </div>
        { modelDropdownOpen.value &&
          <div style={{
            position: 'absolute', top: '100%', left: '0', right: '0', zIndex: '100',
            background: 'var(--white, #fff)', border: '1px solid var(--grey-2)',
            borderRadius: '3px', maxHeight: '300px', overflow: 'auto', boxShadow: '0 2px 8px rgba(0,0,0,0.15)',
          }}>
            <input
              type='text'
              placeholder='Filter models...'
              value={modelSearch.value}
              onInput={(e: Event) => modelSearch.value = (e.target as HTMLInputElement).value}
              style={{width: '100%', padding: '4px 8px', border: 'none', borderBottom: '1px solid var(--grey-2)', outline: 'none'}}
            />
            { filteredModels.value.map((f) =>
              <div
                key={f.id}
                style={{padding: '4px 8px', cursor: 'pointer'}}
                class='hover:bg-slate-100'
                onClick={() => { selectedModel.value = Vue.markRaw(f); modelDropdownOpen.value = false; }}
              >{modelLabel(f)}</div>,
            )}
            { filteredModels.value.length === 0 &&
              <div style={{padding: '4px 8px', color: 'var(--grey-4)'}}>No models found</div> }
          </div>
        }
      </div>
    );

    const renderBasket = () => (
      <div style={{display: 'flex', flexDirection: 'column'}}>
        <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Comparison set</div>
        { entries.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>Select runs in the history and add them here</div> }
        <div class='c2-comparison-rows'>
          { entries.value.map((entry) => {
            const badge = SOURCE_BADGES[entry.sourceKind];
            const isBaseline = entry.id === baselineId.value;
            return <div key={entry.id} class='c2-comparison-row'
              style={{display: 'flex', alignItems: 'center', gap: '6px', padding: '2px 6px'}}>
            <IconFA
              name='thumbtack'
              faStyle={isBaseline ? 'fas' : 'fal'}
              tooltip={isBaseline ? 'Baseline run' : 'Set as baseline'}
              onClick={() => baselineId.value = entry.id}
            />
            <span style={{
              fontSize: '10px', color: 'white', background: badge.color,
              borderRadius: '3px', padding: '0px 4px', flexShrink: '0',
            }}>{badge.label}</span>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', flex: '1'}}
              title={entry.modelName ? `${entry.modelName}: ${entry.name}` : entry.name}>
              {entry.name}
            </span>
            <IconFA name='times' tooltip='Remove from comparison' onClick={() => removeEntry(entry.id)}/>
          </div>;
          })}
        </div>
      </div>
    );

    const renderIndexPickers = () => (
      indexTables.value.length > 0 &&
      <div>
        <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
          <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Index columns</div>
          <div style={{display: 'flex', gap: '12px', alignItems: 'center'}}>
            { renderListFilter(indexFilter, 'Filter tables...') }
            <ToggleInput
              caption='Datetime indexes'
              value={allowDatetimeIndex.value}
              onUpdate:value={(val) => allowDatetimeIndex.value = val}
            />
            <ToggleInput
              caption='Float indexes'
              value={allowFloatIndex.value}
              onUpdate:value={(val) => allowFloatIndex.value = val}
            />
          </div>
        </div>
        <div style={{color: 'var(--grey-4)', paddingBottom: '4px'}}>
          Pick the index (x / key) column for each table to enable column comparison
        </div>
        { indexTables.value.length > 0 && filteredIndexTables.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>No tables match the filter</div> }
        <div class='c2-comparison-rows c2-comparison-table'
          style={{gridTemplateColumns: 'fit-content(480px) max-content 1fr'}}>
          { filteredIndexTables.value.map(({entryId, entryName, table, candidates, current}) => {
            const suggestion = suggestedIndex(candidates);
            return <div key={`${entryId}:${table.path}`} class='c2-comparison-row'
              style={{padding: '2px 6px'}}>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
              title={`${entryName} · ${table.path}`}>
              {entryName} · <span style={{color: 'var(--grey-5)'}}>{table.friendlyPath ?? table.path}</span>
            </span>
            <select
              value={current}
              onChange={(e: Event) => setIndexColumn(entryId, table.path, (e.target as HTMLSelectElement).value)}
              style={{border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 4px', maxWidth: '160px'}}
            >
              <option value=''>— index —</option>
              { candidates.map((col) =>
                <option key={col.name} value={col.name}>
                  {col.name}{col.name === suggestion ? ' (suggested)' : ''}
                </option>,
              )}
            </select>
          </div>;
          })}
        </div>
      </div>
    );

    const renderTargets = () => (
      <div>
        <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
          <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Compare</div>
          { renderListFilter(targetFilter, 'Filter values...') }
        </div>
        { targets.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>
            No comparable data found across the selected runs.
            Add more runs, or set index columns to compare table columns.
          </div> }
        { targets.value.length > 0 && filteredTargets.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>No values match the filter</div> }
        <div class='c2-comparison-rows c2-comparison-table'
          style={{gridTemplateColumns: 'max-content fit-content(400px) max-content max-content max-content 1fr'}}>
          { filteredTargets.value.map((target) => {
            const isSelected = target.key === selectedTargetKey.value;
            return <div
              key={target.key}
              class={isSelected ? 'c2-comparison-row c2-comparison-row-selected' : 'c2-comparison-row'}
              style={{
                padding: '3px 6px', cursor: 'pointer',
                border: `1px solid ${isSelected ? 'var(--blue-1, #2083d5)' : 'transparent'}`,
              }}
              onClick={() => selectedTargetKey.value = target.key}
            >
            <IconFA name={target.kind === 'scalar' ? 'hashtag' : 'table'} tooltip={target.kind}/>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
              title={target.displayName}>
              {target.displayName}
            </span>
            { target.unitsWarning ?
              <IconFA name='exclamation-triangle' tooltip='Units are missing on some runs'/> :
              <span></span> }
            <span style={{
              fontSize: '10px', color: 'white', borderRadius: '3px', padding: '0px 4px',
              background: CONFIDENCE_COLORS[target.confidence], flexShrink: '0',
            }}>{target.confidence}</span>
            <span style={{fontSize: '11px', color: 'var(--grey-4)', flexShrink: '0'}}>
              {target.coverage}/{target.total}
            </span>
          </div>;
          })}
        </div>
      </div>
    );

    const renderStatuses = () => {
      const target = selectedTarget.value;
      const currentComparison = comparison.value;
      if (!target)
        return null;
      const excludedByBuild = new Map(
        (currentComparison?.excluded ?? []).map((item) => [item.entryId, item.reason]));
      return <div style={{display: 'flex', gap: '6px', flexWrap: 'wrap', padding: '4px 0px'}}>
        { entryStatuses.value.map((status) => {
          const entry = entries.value.find((item) => item.id === status.entryId);
          if (!entry)
            return null;
          const buildReason = excludedByBuild.get(status.entryId);
          const matched = status.matched && !buildReason;
          const reason = buildReason ?? status.reason;
          const binding = matched ?
            (target.bindings as {
              entryId: string, path?: string, friendlyPath?: string,
              tablePath?: string, tableFriendlyPath?: string, columnName?: string,
            }[]).find((b) => b.entryId === status.entryId) : undefined;
          const pathText = binding ?
            (binding.friendlyPath ?? binding.path ??
              `${binding.tableFriendlyPath ?? binding.tablePath} · ${binding.columnName}`) : undefined;
          return <span
            key={status.entryId}
            title={pathText ? `${entry.name} · ${pathText}` : entry.name}
            style={{
              fontSize: '11px', borderRadius: '3px', padding: '1px 6px',
              background: matched ? '#e4f3ea' : '#fbeaea',
              color: matched ? '#2c7a4b' : '#a94442',
            }}>
            {entry.name}: {matched ? `matched (${pathText})` : reason}
          </span>;
        })}
      </div>;
    };

    const renderModeToggle = () => (
      <div style={{display: 'flex', gap: '2px', padding: '4px 0px'}}>
        { ([['values', 'Values'], ['delta', 'Δ'], ['deltaPct', 'Δ%']] as [ComparisonMode, string][])
          .map(([value, label]) =>
            <button
              key={value}
              onClick={() => mode.value = value}
              style={{
                padding: '2px 10px', cursor: 'pointer', borderRadius: '3px',
                border: '1px solid var(--grey-2)',
                background: mode.value === value ? 'var(--blue-1, #2083d5)' : 'transparent',
                color: mode.value === value ? 'white' : 'inherit',
              }}
            >{label}</button>,
          )}
      </div>
    );

    const renderComparison = () => {
      const currentComparison = comparison.value;
      if (!selectedTarget.value)
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>Select what to compare above</div>;
      if (!currentComparison || currentComparison.gridDf.rowCount === 0) {
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>
          Nothing to show: no matched data points across the selected runs
        </div>;
      }
      let chart;
      if (currentComparison.kind === 'scalar') {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.gridDf}
          style={{width: '100%', flex: '1', minHeight: '250px'}}
          options={{
            valueColumnName: scalarValueColumn.value,
            valueAggrType: 'avg',
            splitColumnName: RUN_COLUMN,
          }}
        />;
      } else if (currentComparison.isKeyIndex) {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.chartDf}
          style={{width: '100%', flex: '1', minHeight: '250px'}}
          options={{
            valueColumnName: currentComparison.target.displayName,
            valueAggrType: 'avg',
            splitColumnName: currentComparison.indexColumnName,
            stackColumnName: RUN_COLUMN,
          }}
        />;
      } else {
        chart = <Viewer
          type={DG.VIEWER.LINE_CHART}
          dataFrame={currentComparison.chartDf}
          style={{width: '100%', flex: '1', minHeight: '250px'}}
          options={{
            xColumnName: currentComparison.indexColumnName,
            yColumnNames: [currentComparison.target.displayName],
            splitColumnNames: [RUN_COLUMN],
          }}
        />;
      }
      return <div style={{display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px'}}>
        { renderModeToggle() }
        { chart }
        <div style={{height: '220px', paddingTop: '6px'}}>
          <Viewer
            type='Grid'
            dataFrame={currentComparison.gridDf}
            style={{height: '100%', width: '100%'}}
            options={{'allowEdit': false, 'showRowHeader': false}}
          />
        </div>
      </div>;
    };

    return () => Vue.withDirectives(
      <div style={{display: 'flex', height: '100%', width: '100%', overflow: 'hidden'}}>
        <div style={{
          width: `${sidebarWidth.value}px`, flexShrink: '0', display: 'flex', flexDirection: 'column',
          gap: '8px', padding: '8px', overflow: 'auto',
        }}>
          <div style={{fontWeight: 'bold'}}>Model</div>
          { renderModelSelector() }
          { selectedModel.value &&
            <div style={{
              height: `${historyHeight.value}px`, flexShrink: '0',
              border: '1px solid var(--grey-1)', borderRadius: '3px',
            }}>
              <History
                func={selectedModel.value}
                allowCompare={false}
                showIsComplete={true}
                fallbackText='No runs found for this model'
                onSelectionChanged={(calls: DG.FuncCall[]) => historySelection.value = calls.map((c) => Vue.markRaw(c))}
              />
            </div>
          }
          { selectedModel.value &&
            <ResizeHandle
              axis='y'
              size={historyHeight.value}
              min={150}
              max={800}
              onUpdate:size={(size) => historyHeight.value = size}
            />
          }
          { selectedModel.value &&
            <button
              onClick={addSelectedRuns}
              disabled={historySelection.value.length === 0 || isAddingRuns.value}
              style={{
                padding: '4px 8px', cursor: 'pointer', borderRadius: '3px',
                border: '1px solid var(--grey-2)',
              }}
            >
              {isAddingRuns.value ? 'Adding...' :
                `Add selected runs${historySelection.value.length ? ` (${historySelection.value.length})` : ''}`}
            </button>
          }
          <div ref={(el) => {
            if (el && !(el as HTMLElement).contains(tableInput.root))
              (el as HTMLElement).appendChild(tableInput.root);
          }}></div>
          { renderBasket() }
        </div>
        <ResizeHandle
          axis='x'
          size={sidebarWidth.value}
          min={240}
          max={900}
          onUpdate:size={(size) => sidebarWidth.value = size}
        />
        <div style={{flex: '1', display: 'flex', flexDirection: 'column', padding: '8px', overflow: 'auto', minWidth: '0px'}}>
          { entries.value.length < 2 ?
            <div style={{color: 'var(--grey-4)'}}>Add at least two runs (or tables) to compare</div> :
            <div style={{display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px'}}>
              { renderTargets() }
              { renderIndexPickers() }
              { renderStatuses() }
              { renderComparison() }
            </div>
          }
        </div>
      </div>,
      [[ifOverlapping, isLoadingModels.value, 'Loading models...']],
    );
  },
});
