import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {
  IconFA, Viewer, ResizeHandle, ToggleInput, ifOverlapping, wheelGuard,
} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';

import {History} from '../History/History';
import {
  ComparisonTarget, ColumnTarget, MatchConfidence, matchScalarTargets, matchColumnTargets, getEntryStatuses,
  matchesFilter, compatibleTargetsFor, isSplitCandidate, selectionToMap, computeIndexRows,
} from './comparison-core';
import {
  ComparisonEntry, EntrySourceKind, RUN_COLUMN,
  entryFromFuncCall, entryFromDataFrame,
  buildScalarComparison, buildColumnComparison, buildMultiColumnComparison,
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
    const mergeSameFuncs = Vue.ref(true);

    const historySelection = Vue.shallowRef<DG.FuncCall[]>([]);
    const entries = Vue.shallowRef<ComparisonEntry[]>([]);
    // entryId -> tablePath -> index column name
    const indexSelection = Vue.ref<Record<string, Record<string, string>>>({});
    // entryId -> tablePath -> split (category) column name
    const splitSelection = Vue.ref<Record<string, Record<string, string>>>({});
    const selectedTargetKey = Vue.ref('');

    const historyHeight = Vue.ref(320);
    const sidebarWidth = Vue.ref(360);
    // user override for the results chart height; the per-type minimum still applies
    const chartHeight = Vue.ref(0);

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
      // resetting synchronously re-enters the Dart stream controller mid-dispatch
      setTimeout(() => tableInput.value = null);
    });

    const addEntry = (entry: ComparisonEntry) => {
      if (entries.value.some((existing) => existing.id === entry.id))
        return;
      entries.value = [...entries.value, Vue.markRaw(entry)];
    };

    const removeEntry = (id: string) => {
      entries.value = entries.value.filter((entry) => entry.id !== id);
      const {[id]: _removedIndex, ...restIndex} = indexSelection.value;
      indexSelection.value = restIndex;
      const {[id]: _removedSplit, ...restSplit} = splitSelection.value;
      splitSelection.value = restSplit;
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

    const updatedSelection = (
      selection: Record<string, Record<string, string>>,
      members: {entryId: string, tablePath: string}[],
      columnName: string,
    ) => {
      const next = {...selection};
      for (const {entryId, tablePath} of members)
        next[entryId] = {...next[entryId], [tablePath]: columnName};
      return next;
    };

    const setIndexColumn = (members: {entryId: string, tablePath: string}[], columnName: string) => {
      indexSelection.value = updatedSelection(indexSelection.value, members, columnName);
    };

    const setSplitColumn = (members: {entryId: string, tablePath: string}[], columnName: string) => {
      splitSelection.value = updatedSelection(splitSelection.value, members, columnName);
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

    const indexColumnsMap = Vue.computed(() => selectionToMap(indexSelection.value,
      (entryId, tablePath, columnName) => isAllowedIndexType(indexColumnType(entryId, tablePath, columnName))));

    const splitColumnsMap = Vue.computed(() => selectionToMap(splitSelection.value,
      (entryId, tablePath, columnName) => isSplitCandidate(
        {name: columnName, type: indexColumnType(entryId, tablePath, columnName) ?? ''},
        indexColumnsMap.value.get(entryId)?.get(tablePath) ?? '')));

    const targets = Vue.computed<ComparisonTarget[]>(() => {
      const nodes = entries.value.map((entry) => entry.nodes);
      return [
        ...matchScalarTargets(nodes),
        ...matchColumnTargets(nodes, indexColumnsMap.value, splitColumnsMap.value),
      ].sort((a, b) => b.coverage - a.coverage);
    });

    const indexFilter = Vue.ref('');
    const targetFilter = Vue.ref('');

    const filteredTargets = Vue.computed(() => {
      const listed: ComparisonTarget[] = multiMode.value ? compatibleTargets.value : targets.value;
      return listed.filter((target) => matchesFilter(targetFilter.value, target.displayName));
    });

    const selectedTarget = Vue.computed(() =>
      targets.value.find((target) => target.key === selectedTargetKey.value) ?? null);

    const multiMode = Vue.ref(false);
    const multiKeys = Vue.ref<string[]>([]);

    const compatibleTargets = Vue.computed<ColumnTarget[]>(() =>
      compatibleTargetsFor(selectedTarget.value, targets.value, indexColumnType));

    Vue.watch(compatibleTargets, (list) => {
      if (multiMode.value && list.length <= 1)
        multiMode.value = false;
    });

    const chartViewer = Vue.shallowRef<DG.Viewer | null>(null);
    const onChartViewerChanged = (viewer: DG.Viewer | undefined) => {
      chartViewer.value = viewer ? Vue.markRaw(viewer) : null;
    };

    // snapshot export: clone of the chart data plus the chart with its current options
    const openInWorkspace = () => {
      const currentComparison = comparison.value;
      if (!currentComparison || currentComparison.chartDf.rowCount === 0)
        return;
      const nameInput = ui.input.string('Name', {value: currentComparison.chartDf.name || 'Comparison'});
      ui.dialog('Open in workspace')
        .add(nameInput.root)
        .onOK(() => {
          const df = currentComparison.chartDf.clone();
          df.name = nameInput.value || 'Comparison';
          const view = grok.shell.addTableView(df);
          const options = chartViewer.value?.getOptions();
          if (options)
            view.addViewer(options.type, options.look);
        })
        .show({center: true});
    };

    const setMultiMode = (val: boolean) => {
      if (val)
        multiKeys.value = selectedTargetKey.value ? [selectedTargetKey.value] : [];
      else if (multiKeys.value.length > 0)
        selectedTargetKey.value = multiKeys.value[0];
      multiMode.value = val;
    };

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
        const result = buildScalarComparison(target, entries.value);
        return Vue.markRaw({kind: 'scalar' as const, target, ...result});
      }
      if (multiMode.value) {
        const selected = compatibleTargets.value.filter((item) => multiKeys.value.includes(item.key));
        if (selected.length === 0)
          return null;
        if (selected.length > 1) {
          const result = buildMultiColumnComparison(selected, entries.value);
          return result ? Vue.markRaw({kind: 'column' as const, target, ...result}) : null;
        }
        const result = buildColumnComparison(selected[0], entries.value);
        return result ? Vue.markRaw({kind: 'column' as const, target: selected[0], ...result}) : null;
      }
      const result = buildColumnComparison(target, entries.value);
      return result ? Vue.markRaw({kind: 'column' as const, target, ...result}) : null;
    });

    // tables that could participate in column comparison and need an index choice;
    // with merging on, same-function outputs (by nqName) collapse into one row
    const indexRows = Vue.computed(() => computeIndexRows(
      entries.value.map((entry) => entry.nodes),
      indexSelection.value,
      splitSelection.value,
      mergeSameFuncs.value,
      isAllowedIndexType,
    ));

    const filteredIndexRows = Vue.computed(() => indexRows.value.filter((row) =>
      matchesFilter(indexFilter.value, row.entryName ? `${row.entryName} · ${row.label}` : row.label)));

    const suggestedIndex = (columns: {name: string, type: string}[]) =>
      columns.find((col) => col.type === DG.COLUMN_TYPE.DATE_TIME)?.name ??
      columns.find((col) => col.type === DG.COLUMN_TYPE.FLOAT || col.type === DG.COLUMN_TYPE.INT ||
        col.type === DG.COLUMN_TYPE.BIG_INT)?.name;

    const renderColumnSelect = (
      current: string,
      candidates: {name: string}[],
      placeholder: string,
      onSelect: (columnName: string) => void,
      suggestion?: string,
    ) => (
      <select
        value={current}
        onChange={(e: Event) => onSelect((e.target as HTMLSelectElement).value)}
        style={{border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 4px', maxWidth: '160px'}}
      >
        <option value=''>{placeholder}</option>
        { candidates.map((col) =>
          <option key={col.name} value={col.name}>
            {col.name}{col.name === suggestion ? ' (suggested)' : ''}
          </option>,
        )}
      </select>
    );

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
            return <div key={entry.id} class='c2-comparison-row'
              style={{display: 'flex', alignItems: 'center', gap: '6px', padding: '2px 6px'}}>
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
      indexRows.value.length > 0 &&
      <div>
        <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 50px'}}>
          <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Index columns</div>
          <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 12px'}}>
            { renderListFilter(indexFilter, 'Filter tables...') }
            <ToggleInput
              caption='Merge same functions'
              value={mergeSameFuncs.value}
              onUpdate:value={(val) => mergeSameFuncs.value = val}
            />
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
        { filteredIndexRows.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>No tables match the filter</div> }
        <div class='c2-comparison-rows c2-comparison-table'
          style={{gridTemplateColumns: 'fit-content(480px) max-content max-content 1fr'}}>
          { filteredIndexRows.value.map((row) => {
            const suggestion = suggestedIndex(row.candidates);
            return <div key={row.key} class='c2-comparison-row' style={{padding: '2px 6px'}}>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
              title={row.title}>
              { row.coverage ?
                <span style={{fontSize: '11px', color: 'var(--grey-4)'}}>
                  {row.coverage.count}/{row.coverage.total}
                </span> :
                row.entryName
              } · <span style={{color: 'var(--grey-5)'}}>{row.label}</span>
            </span>
            { renderColumnSelect(row.current, row.candidates, '— index —',
              (value) => setIndexColumn(row.members, value), suggestion) }
            { row.splitCandidates.length > 0 ?
              renderColumnSelect(row.currentSplit, row.splitCandidates, '— split —',
                (value) => setSplitColumn(row.members, value)) :
              <span></span> }
          </div>;
          })}
        </div>
      </div>
    );

    const renderTargets = () => (
      <div>
        <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 50px'}}>
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
            const isSelected = multiMode.value ?
              multiKeys.value.includes(target.key) : target.key === selectedTargetKey.value;
            const onRowClick = (e: MouseEvent) => {
              if (!multiMode.value) {
                // grid-like shift+click: select and enter multi-value mode when applicable.
                // Compatible with the current selection -> keep it and add the clicked one;
                // otherwise the clicked target becomes the selection (and the anchor).
                if (e.shiftKey) {
                  if (selectedTargetKey.value &&
                    compatibleTargets.value.some((item) => item.key === target.key)) {
                    multiKeys.value = [...new Set([selectedTargetKey.value, target.key])];
                    multiMode.value = true;
                    return;
                  }
                  selectedTargetKey.value = target.key;
                  if (compatibleTargetsFor(target, targets.value, indexColumnType).length > 1) {
                    multiKeys.value = [target.key];
                    multiMode.value = true;
                  }
                  return;
                }
                selectedTargetKey.value = target.key;
                return;
              }
              multiKeys.value = multiKeys.value.includes(target.key) ?
                multiKeys.value.filter((key) => key !== target.key) : [...multiKeys.value, target.key];
            };
            return <div
              key={target.key}
              class={isSelected ? 'c2-comparison-row c2-comparison-row-selected' : 'c2-comparison-row'}
              style={{
                padding: '3px 6px', cursor: 'pointer',
                borderLeft: `3px solid ${isSelected ? 'var(--blue-1, #2083d5)' : 'transparent'}`,
              }}
              onClick={onRowClick}
              onMousedown={(e: MouseEvent) => { if (e.shiftKey) e.preventDefault(); }}
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

    const renderResultsHeader = () => (
      <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 50px'}}>
        <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Results</div>
        <div style={{display: 'flex', alignItems: 'center', gap: '4px 12px', flexWrap: 'wrap'}}>
          { compatibleTargets.value.length > 1 &&
            <ToggleInput
              caption='Multiple values'
              value={multiMode.value}
              onUpdate:value={setMultiMode}
            />
          }
          { comparison.value && comparison.value.chartDf.rowCount > 0 &&
            <span
              style={{
                display: 'flex', alignItems: 'center', gap: '4px',
                cursor: 'pointer', color: 'var(--blue-1, #2083d5)',
              }}
              onClick={openInWorkspace}
            >
              <IconFA
                name='external-link'
                tooltip='Open a snapshot of the data and chart in the workspace'
              />
              <span>Open in workspace</span>
            </span>
          }
        </div>
      </div>
    );

    const renderStatuses = () => {
      const target = selectedTarget.value;
      if (!target)
        return null;
      // only runs that could not participate are flagged; matched runs are visible in the chart
      const problems = entryStatuses.value
        .map((status) => {
          const entry = entries.value.find((item) => item.id === status.entryId);
          if (!entry || status.matched)
            return null;
          return {entry, reason: status.reason};
        })
        .filter((item) => item != null);
      if (problems.length === 0)
        return null;
      return <div style={{display: 'flex', gap: '6px', flexWrap: 'wrap', padding: '4px 0px'}}>
        { problems.map(({entry, reason}) =>
          <span
            key={entry.id}
            title={entry.name}
            style={{
              fontSize: '11px', borderRadius: '3px', padding: '1px 6px',
              background: '#fbeaea', color: '#a94442',
            }}>
            {entry.name}: {reason}
          </span>,
        )}
      </div>;
    };

    const renderComparison = () => {
      const currentComparison = comparison.value;
      if (!selectedTarget.value)
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>Select what to compare above</div>;
      if (!currentComparison || currentComparison.chartDf.rowCount === 0) {
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>
          Nothing to show: no data points across the selected runs
        </div>;
      }
      const chartMinHeight = 'valueColumnNames' in currentComparison ?
        250 * Math.max(1, (currentComparison.valueColumnNames as string[]).length) : 250;
      const effectiveChartHeight = Math.max(chartMinHeight, chartHeight.value);
      const chartStyle = {width: '100%', height: `${effectiveChartHeight}px`, flexShrink: '0'};
      let chart;
      if (currentComparison.kind === 'scalar') {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={{
            valueColumnName: currentComparison.target.displayName,
            valueAggrType: 'avg',
            splitColumnName: RUN_COLUMN,
          }}
        />;
      } else if (currentComparison.isKeyIndex && !('valueColumnNames' in currentComparison)) {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={{
            valueColumnName: currentComparison.target.displayName,
            valueAggrType: 'avg',
            splitColumnName: currentComparison.indexColumnName,
            stackColumnName: RUN_COLUMN,
          }}
        />;
      } else {
        const yColumnNames = 'valueColumnNames' in currentComparison ?
          currentComparison.valueColumnNames as string[] : [currentComparison.target.displayName];
        chart = <Viewer
          type={DG.VIEWER.LINE_CHART}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={{
            xColumnName: currentComparison.indexColumnName,
            yColumnNames,
            splitColumnNames: currentComparison.splitColumnName ?
              [RUN_COLUMN, currentComparison.splitColumnName] : [RUN_COLUMN],
            multiAxis: false,
          }}
        />;
      }
      return <div style={{display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px'}}>
        <ResizeHandle
          axis='y'
          size={effectiveChartHeight}
          min={chartMinHeight}
          reverse={true}
          onUpdate:size={(size) => chartHeight.value = size}
        />
        { Vue.withDirectives(chart, [[wheelGuard]]) }
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
            <div style={{display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px', gap: '12px'}}>
              { renderIndexPickers() }
              { renderTargets() }
              { renderResultsHeader() }
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
