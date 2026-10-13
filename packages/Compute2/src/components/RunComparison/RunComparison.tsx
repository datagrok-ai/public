import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {onClickOutside} from '@vueuse/core';

import {
  IconFA, Viewer, ResizeHandle, ToggleInput, ifOverlapping, wheelGuard,
} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getRunTitle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';

import {History} from '../History/History';
import {
  ComparisonTarget, ColumnTarget, ColumnCandidate, CandidateOverrides, MatchConfidence,
  ComparisonEntry, EntrySourceKind, NameMapping, UnitMapping, candidateId,
  AxisModeSelection, AxisModeConfig, AxisMode, TimeUnit, TIME_UNITS, isIndexCandidateType, isNumericType,
} from './types';
import {
  matchScalarTargets, matchColumnTargets, buildAliasGroups, buildUnitConverter,
} from './matching';
import {
  getEntryStatuses, matchesFilter, compatibleTargetsFor, selectionGaps,
  isSplitCandidate, selectionToMap, computeIndexRows, isTargetEqualAcrossRuns, resolveAxisModes,
  targetContexts, EntryStatus,
} from './selection';
import {entryFromFuncCall, entryFromDataFrame} from './entry-extraction';
import {
  TargetSelection, EMPTY_TARGET_SELECTION, isCharted, clickTarget, setMultiMode, reconcileTargets, pruneOverrides,
} from './target-selection';
import {canSaveProject, captureTableLayout, showProjectSaveDialog} from '../../project-export';
import {buildComparison, columnPrimitives} from './comparison-builders';
import {showMappingsDialog} from './mappings-dialog';
import {ScalarChartType, chartSpec} from './chart-spec';
import './RunComparison.css';

// chip styling lives in RunComparison.css (c2-comparison-badge + these modifiers)
const sourceBadge = (kind: EntrySourceKind) =>
  <span class={`c2-comparison-badge c2-badge-${kind}`}>{kind}</span>;

const WARNING_MODES: Record<NonNullable<EntryStatus['warning']>, string> = {
  'relative timeseries not set': 'relative timeseries',
  'independent points not set': 'independent points',
};

const confidenceBadge = (confidence: MatchConfidence) =>
  <span class={`c2-comparison-badge c2-confidence-${confidence}`}>{confidence}</span>;

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

    const mergeSameFuncs = Vue.ref(true);
    const hideEqual = Vue.ref(true);

    const historySelection = Vue.shallowRef<DG.FuncCall[]>([]);
    const entries = Vue.shallowRef<ComparisonEntry[]>([]);
    // entryId -> tablePath -> index column name
    const indexSelection = Vue.ref<Record<string, Record<string, string>>>({});
    // entryId -> tablePath -> split (category) column name
    const splitSelection = Vue.ref<Record<string, Record<string, string>>>({});
    // entryId -> tablePath -> axis-mode pick; stale picks are ignored, not deleted
    const axisModeSelection = Vue.ref<AxisModeSelection>({});
    const targetSelection = Vue.shallowRef<TargetSelection>(EMPTY_TARGET_SELECTION);
    const selectedTargetKey = Vue.computed(() => targetSelection.value.selectedKey);
    const multiMode = Vue.computed(() => targetSelection.value.multiMode);
    const multiKeys = Vue.computed(() => targetSelection.value.multiKeys);
    const candidateOverrides = Vue.ref<CandidateOverrides>({});
    const expandedTargetKeys = Vue.ref<Record<string, boolean>>({});

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

    onClickOutside(modelDropdownRef, () => modelDropdownOpen.value = false);

    const modelLabel = (f: DG.Func) => f.friendlyName ?? f.name;
    const filteredModels = Vue.computed(() => {
      const query = modelSearch.value.toLowerCase();
      return models.value.filter((f) => !query || modelLabel(f).toLowerCase().includes(query));
    });

    const tableInput = ui.input.table('Add table');
    const tableInputSub = tableInput.onChanged.subscribe(() => {
      const df = tableInput.value;
      if (!df)
        return;
      addEntry(entryFromDataFrame(df));
      // resetting synchronously re-enters the Dart stream controller mid-dispatch
      setTimeout(() => tableInput.value = null);
    });

    Vue.onBeforeUnmount(() => {
      tableInputSub.unsubscribe();
      // input.root never enters this component's DOM tree, so nothing else disposes it
      ui.remove(tableInput.root);
    });

    // pre-fill pickers from {comparisonIndex}/{comparisonSplit} output annotations;
    // stored (user) selections always win
    const applyAnnotatedDefaults = (entry: ComparisonEntry) => {
      for (const table of entry.nodes.tables) {
        const member = [{entryId: entry.id, tablePath: table.path}];
        const indexDefault = table.defaultIndexColumn;
        if (indexDefault && !indexSelection.value[entry.id]?.[table.path] &&
          table.columns.some((col) => col.name === indexDefault))
          setIndexColumn(member, indexDefault);
        const splitDefault = table.defaultSplitColumn;
        const splitColumn = table.columns.find((col) => col.name === splitDefault);
        if (splitDefault && !splitSelection.value[entry.id]?.[table.path] &&
          splitColumn && isSplitCandidate(splitColumn, indexSelection.value[entry.id]?.[table.path] ?? ''))
          setSplitColumn(member, splitDefault);
        const modeDefault = table.defaultAxisMode;
        if (modeDefault && !axisModeSelection.value[entry.id]?.[table.path]) {
          setAxisMode(member, {
            mode: modeDefault,
            // units matter only on the elapsed axis; storing them elsewhere would let a
            // unit the user never saw outrank the column-meta prefill later
            ...modeDefault === 'timeseries' && table.defaultTimeUnits ?
              {units: table.defaultTimeUnits} : {},
          });
        }
      }
    };

    const addEntry = (entry: ComparisonEntry) => {
      if (entries.value.some((existing) => existing.id === entry.id))
        return;
      entries.value = [...entries.value, Vue.markRaw(entry)];
      applyAnnotatedDefaults(entry);
    };

    const removeEntry = (id: string) => {
      entries.value = entries.value.filter((entry) => entry.id !== id);
      const {[id]: _removedIndex, ...restIndex} = indexSelection.value;
      indexSelection.value = restIndex;
      const {[id]: _removedSplit, ...restSplit} = splitSelection.value;
      splitSelection.value = restSplit;
      const {[id]: _removedAxis, ...restAxis} = axisModeSelection.value;
      axisModeSelection.value = restAxis;
    };

    const addSelectedRuns = async () => {
      if (isAddingRuns.value || historySelection.value.length === 0)
        return;
      isAddingRuns.value = true;
      try {
        const calls = historySelection.value.filter((call) => !entriesById.value.has(call.id));
        const results = await Promise.allSettled(calls.map(async (call) =>
          entryFromFuncCall(await historyUtils.loadRun(call.id))));
        const failedRuns: string[] = [];
        results.forEach((result, i) => {
          if (result.status === 'rejected') {
            console.error(`Run comparison: failed to load run ${calls[i].id}`, result.reason);
            failedRuns.push(getRunTitle(calls[i]));
            return;
          }
          addEntry(result.value);
          if (result.value.failedSteps)
            grok.shell.warning(`${result.value.name}: could not load ${result.value.failedSteps.join(', ')}`);
        });
        if (failedRuns.length > 0)
          grok.shell.error(`Could not load ${failedRuns.join(', ')}`);
      } finally {
        isAddingRuns.value = false;
      }
    };

    const updatedSelection = <T, >(
      selection: Record<string, Record<string, T>>,
      members: {entryId: string, tablePath: string}[],
      value: T,
    ) => {
      const next = {...selection};
      for (const {entryId, tablePath} of members)
        next[entryId] = {...next[entryId], [tablePath]: value};
      return next;
    };

    const setIndexColumn = (members: {entryId: string, tablePath: string}[], columnName: string) => {
      indexSelection.value = updatedSelection(indexSelection.value, members, columnName);
    };

    const setSplitColumn = (members: {entryId: string, tablePath: string}[], columnName: string) => {
      splitSelection.value = updatedSelection(splitSelection.value, members, columnName);
    };

    const setAxisMode = (members: {entryId: string, tablePath: string}[], config: AxisModeConfig) => {
      axisModeSelection.value = updatedSelection(axisModeSelection.value, members, config);
    };

    const entriesById = Vue.computed(() => new Map(entries.value.map((entry) => [entry.id, entry])));

    const tableNode = (entryId: string, tablePath: string) =>
      entriesById.value.get(entryId)?.nodes.tables.find((table) => table.path === tablePath);

    const indexColumnType = (entryId: string, tablePath: string, columnName: string) =>
      tableNode(entryId, tablePath)?.columns.find((col) => col.name === columnName)?.type;

    const defaultIndexOf = (entryId: string, tablePath: string) => tableNode(entryId, tablePath)?.defaultIndexColumn;

    const indexColumnsMap = Vue.computed(() => selectionToMap(indexSelection.value,
      (entryId, tablePath, columnName) =>
        isIndexCandidateType(indexColumnType(entryId, tablePath, columnName)) ||
        columnName === defaultIndexOf(entryId, tablePath)));

    const splitColumnsMap = Vue.computed(() => selectionToMap(splitSelection.value,
      (entryId, tablePath, columnName) => isSplitCandidate(
        {name: columnName, type: indexColumnType(entryId, tablePath, columnName) ?? ''},
        indexColumnsMap.value.get(entryId)?.get(tablePath) ?? '')));

    const nameMappings = Vue.ref<NameMapping[]>([]);
    const aliasGroups = Vue.computed(() => buildAliasGroups(nameMappings.value));
    const unitMappings = Vue.ref<UnitMapping[]>([]);
    const unitConverter = Vue.computed(() => buildUnitConverter(unitMappings.value));
    const mappingCount = Vue.computed(() => nameMappings.value.length + unitMappings.value.length);

    const targets = Vue.computed<ComparisonTarget[]>(() => {
      const nodes = entries.value.map((entry) => entry.nodes);
      return [
        ...matchScalarTargets(nodes, aliasGroups.value, unitConverter.value),
        ...matchColumnTargets(nodes, indexColumnsMap.value, splitColumnsMap.value, candidateOverrides.value,
          aliasGroups.value, unitConverter.value),
      ].sort((a, b) => a.kind !== b.kind ?
        (a.kind === 'scalar' ? -1 : 1) :
        a.displayName.localeCompare(b.displayName));
    });

    const contexts = Vue.computed(() => targetContexts(targets.value));

    // radio semantics: checking an item makes it the run's pick and unchecks the
    // run's current one; unchecking excludes the run from the target
    const setCandidateEnabled = (target: ColumnTarget, candidate: ColumnCandidate, enabled: boolean) => {
      const next = {...candidateOverrides.value[target.key]};
      if (enabled) {
        for (const other of target.candidates) {
          if (other !== candidate && other.enabled && other.binding.entryId === candidate.binding.entryId)
            next[candidateId(other.binding)] = false;
        }
      }
      next[candidateId(candidate.binding)] = enabled;
      candidateOverrides.value = {...candidateOverrides.value, [target.key]: next};
    };

    const resetTargetOverrides = (targetKey: string) => {
      const {[targetKey]: _removed, ...rest} = candidateOverrides.value;
      candidateOverrides.value = rest;
    };

    const toggleExpanded = (targetKey: string) => {
      expandedTargetKeys.value = {
        ...expandedTargetKeys.value,
        [targetKey]: !expandedTargetKeys.value[targetKey],
      };
    };

    // removed entries are not consumed by matching, so the override rewrite converges in one cycle
    Vue.watch(targets, (list) => {
      candidateOverrides.value = pruneOverrides(candidateOverrides.value, list);
      targetSelection.value = reconcileTargets(targetSelection.value, new Set(list.map((target) => target.key)));
    });

    const indexFilter = Vue.ref('');
    const targetFilter = Vue.ref('');

    const columnValues = (entryId: string, tablePath: string, columnName: string) => {
      const col = entriesById.value.get(entryId)?.dataFrames.get(tablePath)?.col(columnName);
      return col ? columnPrimitives(col) : undefined;
    };

    const equalTargetKeys = Vue.computed(() => {
      if (!hideEqual.value)
        return new Set<string>();
      const cache = new Map<string, unknown[] | undefined>();
      const cached = (entryId: string, tablePath: string, columnName: string) => {
        const key = `${entryId}|${tablePath}|${columnName}`;
        if (!cache.has(key))
          cache.set(key, columnValues(entryId, tablePath, columnName));
        return cache.get(key);
      };
      return new Set(targets.value
        .filter((target) => isTargetEqualAcrossRuns(target, cached, indexColumnType))
        .map((target) => target.key));
    });

    // in multi mode: suggestions plus everything already selected — a selected target
    // that became conflicting must stay visible so it can be unchecked. A charted value
    // stays listed even when equal, so what is on the chart can always be deselected
    const filteredTargets = Vue.computed(() => {
      const charted = (target: ComparisonTarget) => isCharted(targetSelection.value, target.key);
      const listed: ComparisonTarget[] = multiMode.value ?
        targets.value.filter((target) =>
          charted(target) || compatibleTargets.value.some((item) => item.key === target.key)) :
        targets.value;
      return listed.filter((target) =>
        (!equalTargetKeys.value.has(target.key) || charted(target)) &&
        matchesFilter(targetFilter.value, `${target.displayName} ${contexts.value.get(target.key) ?? ''}`));
    });

    const selectedTarget = Vue.computed(() =>
      targets.value.find((target) => target.key === selectedTargetKey.value) ?? null);

    const scalarChartType = Vue.ref<ScalarChartType>('radar');
    const radarAvailable = DG.Func.find({meta: {role: DG.FUNC_TYPES.VIEWER}})
      .some((f) => f.friendlyName === 'Radar');

    const compatibleTargets = Vue.computed(() =>
      compatibleTargetsFor(selectedTarget.value, targets.value, indexColumnType));

    const chartViewer = Vue.shallowRef<DG.Viewer | null>(null);
    const onChartViewerChanged = (viewer: DG.Viewer | undefined) => {
      chartViewer.value = viewer ? Vue.markRaw(viewer) : null;
      // JsViewer roots (ui.box) expect the host to size them; without this the
      // radar collapses to 0x0 and echarts renders an empty chart
      if (viewer) {
        viewer.root.style.width = '100%';
        viewer.root.style.height = '100%';
      }
    };

    // snapshot export: clone of the chart data plus the chart with its current options;
    // null when the comparison went away behind the non-modal export dialog
    const snapshotData = (name: string) => {
      const current = comparison.value;
      if (!current)
        return null;
      const df = current.chartDf.clone();
      df.name = name;
      const options = chartViewer.value?.getOptions();
      return {df, viewers: options ? [{type: options.type as string, options: options.look}] : []};
    };

    const snapshotView = (name: string) => {
      const snap = snapshotData(name);
      if (!snap)
        return;
      const {df, viewers} = snap;
      const view = grok.shell.addTableView(df);
      for (const {type, options} of viewers)
        view.addViewer(type, options);
      return view;
    };

    const saveAndShare = async (name: string) => {
      const snap = snapshotData(name);
      if (!snap)
        return;
      const {df, viewers} = snap;
      const layout = captureTableLayout(df, viewers);
      try {
        await showProjectSaveDialog([df], [layout], name);
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        grok.shell.closeTable(df);
      }
    };

    const exportComparison = () => {
      const currentComparison = comparison.value;
      if (!currentComparison || currentComparison.chartDf.rowCount === 0)
        return;
      const nameInput = ui.input.string('Name', {value: currentComparison.chartDf.name || 'Comparison'});
      const dlg = ui.dialog('Export comparison').add(nameInput.root);
      const openLabel = 'OPEN IN WORKSPACE';
      if (canSaveProject()) {
        dlg.addButton(openLabel, () => {
          dlg.close();
          snapshotView(nameInput.value || 'Comparison');
        });
        dlg.onOK(() => saveAndShare(nameInput.value || 'Comparison'));
        dlg.show({center: true});
        dlg.getButton('OK').innerText = 'SAVE & SHARE';
      } else {
        dlg.onOK(() => snapshotView(nameInput.value || 'Comparison'));
        dlg.show({center: true});
        dlg.getButton('OK').innerText = openLabel;
      }
    };

    const editMappings = () => showMappingsDialog(nameMappings.value, unitMappings.value, (names, units) => {
      nameMappings.value = names;
      unitMappings.value = units;
    });

    // DG inputs fire onChanged on programmatic sets too, so rendering a new mode value
    // echoes back through the toggle — the no-op transition keeps multiKeys intact
    const onMultiModeChanged = (on: boolean) => {
      targetSelection.value = setMultiMode(targetSelection.value, on);
    };

    const compatibleKeys = (anchorKey: string) => compatibleTargetsFor(
      targets.value.find((target) => target.key === anchorKey) ?? null, targets.value, indexColumnType,
    ).map((target) => target.key);

    // non-default axis modes only; tables whose index moved off a time type drop out
    const resolvedAxisModes = Vue.computed(() => resolveAxisModes(
      entries.value.map((entry) => entry.nodes),
      indexColumnsMap.value,
      axisModeSelection.value,
    ));

    const entryStatuses = Vue.computed(() => getEntryStatuses(
      entries.value.map((entry) => entry.nodes),
      selectedTarget.value,
      indexColumnsMap.value,
      {axisModes: resolvedAxisModes.value, aliases: aliasGroups.value, convert: unitConverter.value},
    ));

    // the multi-value selection in chart order: the anchor first, then list order
    const selectedTargets = <T extends ComparisonTarget>(anchor: T): T[] => {
      const selected = targets.value.filter((item): item is T =>
        item.kind === anchor.kind && multiKeys.value.includes(item.key));
      const anchorIndex = selected.findIndex((item) => item.key === anchor.key);
      if (anchorIndex > 0)
        selected.unshift(...selected.splice(anchorIndex, 1));
      return selected;
    };

    const gaps = Vue.computed(() => multiMode.value && selectedTarget.value ?
      selectionGaps(selectedTargets(selectedTarget.value)) : new Map<string, string[]>());

    const comparison = Vue.computed(() => {
      const anchor = selectedTarget.value;
      if (!anchor)
        return null;
      const result = buildComparison(
        multiMode.value ? selectedTargets(anchor) : [anchor], entries.value, resolvedAxisModes.value);
      return result ? Vue.markRaw(result) : null;
    });

    const currentChartSpec = Vue.computed(() =>
      chartSpec(comparison.value, selectedTarget.value != null, scalarChartType.value, radarAvailable));

    // tables that could participate in column comparison and need an index choice;
    // with merging on, same-function outputs (by nqName) collapse into one row
    const indexRows = Vue.computed(() => computeIndexRows(
      entries.value.map((entry) => entry.nodes),
      indexSelection.value,
      splitSelection.value,
      mergeSameFuncs.value,
      axisModeSelection.value,
    ));

    const filteredIndexRows = Vue.computed(() => indexRows.value.filter((row) =>
      matchesFilter(indexFilter.value, row.entryName ? `${row.entryName} · ${row.label}` : row.label)));

    const suggestedIndex = (columns: {name: string, type: string}[]) =>
      columns.find((col) => col.type === DG.COLUMN_TYPE.DATE_TIME)?.name ??
      columns.find((col) => isNumericType(col.type))?.name;

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
            background: 'var(--white)', border: '1px solid var(--grey-2)',
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
          { entries.value.map((entry) =>
            <div key={entry.id} class='c2-comparison-row'
              style={{display: 'flex', alignItems: 'center', gap: '6px', padding: '2px 6px'}}>
              { sourceBadge(entry.sourceKind) }
              <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', flex: '1'}}
                title={entry.modelName ? `${entry.modelName}: ${entry.name}` : entry.name}>
                {entry.name}
              </span>
              <IconFA name='times' tooltip='Remove from comparison' onClick={() => removeEntry(entry.id)}/>
            </div>,
          )}
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
          </div>
        </div>
        <div style={{color: 'var(--grey-4)', paddingBottom: '4px'}}>
          Pick the index (x / key) column for each table to enable column comparison
        </div>
        { filteredIndexRows.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>No tables match the filter</div> }
        <div class='c2-comparison-rows c2-comparison-table'
          style={{gridTemplateColumns: 'fit-content(480px) max-content max-content max-content max-content 1fr'}}>
          { filteredIndexRows.value.map((row) => {
            const suggestion = suggestedIndex(row.candidates);
            const axis = row.axis;
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
            { axis ?
              <select
                value={axis.mode}
                onChange={(e: Event) => setAxisMode(row.members,
                  {mode: (e.target as HTMLSelectElement).value as AxisMode, units: axis.units})}
                title={'How the table rows relate across runs: an indexed series, a timeseries ' +
                  'relative to each run\'s start, or independent points'}
                style={{border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 4px'}}
              >
                <option value='series'>series</option>
                <option value='timeseries'>relative timeseries</option>
                <option value='points'>independent points</option>
              </select> :
              <span></span> }
            { axis?.mode === 'timeseries' && axis.units != null ?
              <select
                value={axis.units}
                onChange={(e: Event) => setAxisMode(row.members,
                  {mode: 'timeseries', units: (e.target as HTMLSelectElement).value as TimeUnit})}
                title='Time units of the index column values'
                style={{border: '1px solid var(--grey-2)', borderRadius: '3px', padding: '1px 4px'}}
              >
                { TIME_UNITS.map((unit) =>
                  <option key={unit} value={unit}>{unit}</option>,
                )}
              </select> :
              <span></span> }
          </div>;
          })}
        </div>
      </div>
    );

    const renderCandidateRow = (target: ColumnTarget, candidate: ColumnCandidate) => {
      const id = candidateId(candidate.binding);
      return <label
        key={id}
        style={{display: 'flex', alignItems: 'center', gap: '6px', padding: '1px 0px 1px 14px', cursor: 'pointer'}}
      >
        <input
          type='checkbox'
          checked={candidate.enabled}
          onChange={(e: Event) => setCandidateEnabled(target, candidate, (e.target as HTMLInputElement).checked)}
        />
        <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
          title={candidate.binding.tablePath}>
          {candidate.binding.tableFriendlyPath ?? candidate.binding.tableName}
        </span>
        <span style={{color: 'var(--grey-5)', flexShrink: '0'}}>{candidate.binding.columnName}</span>
        { confidenceBadge(candidate.confidence) }
        { candidate.unitsWarn &&
          <IconFA name='exclamation-triangle' tooltip='Units are missing on some items'/> }
        { candidate.auto &&
          <span style={{fontSize: '11px', color: 'var(--grey-4)', flexShrink: '0'}}>(auto)</span> }
      </label>;
    };

    const renderCandidatePanel = (target: ColumnTarget) => {
      const groups: {entry: ComparisonEntry, candidates: ColumnCandidate[]}[] = [];
      for (const candidate of target.candidates) {
        const last = groups[groups.length - 1];
        if (last && last.entry.id === candidate.binding.entryId) {
          last.candidates.push(candidate);
          continue;
        }
        const entry = entriesById.value.get(candidate.binding.entryId);
        if (entry)
          groups.push({entry, candidates: [candidate]});
      }
      return <div
        key={`${target.key}:candidates`}
        class='c2-comparison-expansion'
        style={{gridColumn: '1 / -1', padding: '2px 6px 6px 26px'}}
        onClick={(e: MouseEvent) => e.stopPropagation()}
      >
        { candidateOverrides.value[target.key] &&
          <span
            style={{fontSize: '11px', color: 'var(--blue-1)', cursor: 'pointer'}}
            onClick={() => resetTargetOverrides(target.key)}
          >Reset to auto</span> }
        { groups.map(({entry, candidates}) =>
          <div key={entry.id}>
            <div style={{display: 'flex', alignItems: 'center', gap: '6px', padding: '2px 0px'}}>
              { sourceBadge(entry.sourceKind) }
              <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
                title={entry.name}>{entry.name}</span>
            </div>
            { candidates.map((candidate) => renderCandidateRow(target, candidate)) }
          </div>,
        )}
      </div>;
    };

    const renderTargets = () => (
      <div>
        <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 50px'}}>
          <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Compare</div>
          <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 12px'}}>
            { renderListFilter(targetFilter, 'Filter values...') }
            <ToggleInput
              caption='Hide equal values'
              value={hideEqual.value}
              onUpdate:value={(val) => hideEqual.value = val}
            />
            { (multiMode.value || compatibleTargets.value.length > 1) &&
              <ToggleInput
                caption='Multiple values'
                tooltip='Chart several values together; Shift+click on a row does the same'
                value={multiMode.value}
                onUpdate:value={onMultiModeChanged}
              />
            }
            <div style={{display: 'flex', alignItems: 'center', gap: '4px'}}>
              <span class='c2-comparison-link' onClick={editMappings}>
                <IconFA name='exchange-alt' tooltip='Names that mean the same value, and unit conversion factors'/>
                <span>Mappings...</span>
              </span>
              { mappingCount.value > 0 &&
                <span style={{fontSize: '11px', color: 'var(--grey-4)'}}>{mappingCount.value}</span> }
            </div>
          </div>
        </div>
        { targets.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>
            No comparable data found across the selected runs.
            Add more runs, or set index columns to compare table columns.
          </div> }
        { targets.value.length > 0 && filteredTargets.value.length === 0 &&
          <div style={{color: 'var(--grey-4)'}}>
            { !targetFilter.value && equalTargetKeys.value.size > 0 ?
              'All values are equal across the selected runs' : 'No values match the filter' }
          </div> }
        <div class='c2-comparison-rows c2-comparison-table'
          style={{
            gridTemplateColumns: 'max-content max-content fit-content(400px) max-content max-content max-content 1fr',
          }}>
          { filteredTargets.value.map((target) => {
            const isSelected = isCharted(targetSelection.value, target.key);
            const onRowClick = (e: MouseEvent) => {
              targetSelection.value = clickTarget(targetSelection.value, target.key, e.shiftKey, compatibleKeys);
            };
            const isExpanded = target.kind === 'column' && !!expandedTargetKeys.value[target.key];
            const context = contexts.value.get(target.key);
            const gapRuns = (gaps.value.get(target.key) ?? [])
              .map((id) => entriesById.value.get(id)?.name ?? id);
            const row = <div
              key={target.key}
              class={isSelected ? 'c2-comparison-row c2-comparison-row-selected' : 'c2-comparison-row'}
              style={{
                padding: '3px 6px', cursor: 'pointer',
                borderLeft: `3px solid ${isSelected ? 'var(--blue-1)' : 'transparent'}`,
              }}
              onClick={onRowClick}
              onMousedown={(e: MouseEvent) => { if (e.shiftKey) e.preventDefault(); }}
            >
            { target.kind === 'column' ?
              <span
                class='c2-comparison-expand'
                onClick={(e: Event) => {e.stopPropagation(); toggleExpanded(target.key);}}
              >
                <IconFA name={isExpanded ? 'chevron-down' : 'chevron-right'} tooltip='Show matched items'/>
              </span> :
              <span></span> }
            <IconFA name={target.kind === 'scalar' ? 'hashtag' : 'table'} tooltip={target.kind}/>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
              title={context ? `${target.displayName} · ${context}` : target.displayName}>
              {target.displayName}
              { context && <span class='c2-comparison-context'>{context}</span> }
            </span>
            <span class='c2-comparison-units' title={target.units ? 'Units of the first run' : ''}>{target.units ?? ''}</span>
            { target.unitsWarning ?
              <IconFA name='exclamation-triangle' tooltip='Units are missing on some runs'/> :
              <span></span> }
            { confidenceBadge(target.confidence) }
            <span style={{fontSize: '11px', color: 'var(--grey-4)', flexShrink: '0'}}>
              {target.coverage}/{target.total}
            </span>
            { gapRuns.length > 0 ?
              <span
                class='c2-comparison-badge c2-badge-partial'
                title={target.kind === 'scalar' ?
                  `No value on: ${gapRuns.join(', ')}` :
                  `Shown as gaps (no data in the charted rows): ${gapRuns.join(', ')}`}
              >partial</span> :
              <span></span> }
          </div>;
            return isExpanded ? [row, renderCandidatePanel(target as ColumnTarget)] : row;
          })}
        </div>
      </div>
    );

    const renderResultsHeader = () => (
      <div style={{display: 'flex', alignItems: 'center', flexWrap: 'wrap', gap: '4px 50px'}}>
        <div style={{fontWeight: 'bold', padding: '4px 0px'}}>Results</div>
        <div style={{display: 'flex', alignItems: 'center', gap: '4px 12px', flexWrap: 'wrap'}}>
          { comparison.value && comparison.value.chartDf.rowCount > 0 &&
            <span class='c2-comparison-link' onClick={exportComparison}>
              <IconFA
                name='external-link'
                tooltip='Open a snapshot of the data and chart in the workspace, or save and share it as a project'
              />
              <span>Export...</span>
            </span>
          }
        </div>
      </div>
    );

    const renderStatuses = () => {
      const target = selectedTarget.value;
      if (!target)
        return null;
      // excluded runs get red chips; matched runs on a partially-configured axis-mode
      // target get amber ones — they still chart, just as a plain series
      const chips = entryStatuses.value
        .map((status) => {
          const entry = entriesById.value.get(status.entryId);
          if (!entry)
            return null;
          if (!status.matched)
            return {entry, text: status.units ? `${status.reason} (${status.units})` : status.reason};
          if (status.warning)
            return {entry, text: status.warning, warningMode: WARNING_MODES[status.warning]};
          return null;
        })
        .filter((item) => item != null);
      if (chips.length === 0)
        return null;
      return <div style={{display: 'flex', gap: '6px', flexWrap: 'wrap', padding: '4px 0px'}}>
        { chips.map(({entry, text, warningMode}) =>
          <span
            key={entry.id}
            class={`c2-comparison-status ${warningMode ? 'c2-status-warning' : 'c2-status-error'}`}
            title={warningMode ?
              `The value charts as a plain series until every participating table has the ${warningMode} mode set` :
              entry.name}>
            {entry.name}: {text}
          </span>,
        )}
      </div>;
    };

    const renderComparison = () => {
      const spec = currentChartSpec.value;
      if (spec.kind === 'hint')
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>{spec.text}</div>;
      const effectiveChartHeight = Math.max(spec.minHeight, chartHeight.value);
      return <div style={{display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px'}}>
        <ResizeHandle
          axis='y'
          size={effectiveChartHeight}
          min={spec.minHeight}
          reverse={true}
          onUpdate:size={(size) => chartHeight.value = size}
        />
        { spec.scalarChartSwitch &&
          <div style={{display: 'flex', gap: '4px', padding: '2px 0px', flexShrink: '0'}}>
            { ([['radar', 'Radar'], ['pcplot', 'PC Plot']] as const).map(([value, label]) =>
              <span
                key={value}
                onClick={() => scalarChartType.value = value}
                style={{
                  cursor: 'pointer', fontSize: '11px', borderRadius: '3px', padding: '1px 8px',
                  border: '1px solid var(--grey-2)',
                  background: scalarChartType.value === value ? 'var(--blue-1)' : 'transparent',
                  color: scalarChartType.value === value ? 'var(--white)' : 'var(--grey-5)',
                }}
              >{label}</span>,
            )}
          </div> }
        { Vue.withDirectives(
          <Viewer
            type={spec.type}
            dataFrame={comparison.value!.chartDf}
            style={{width: '100%', height: `${effectiveChartHeight}px`, flexShrink: '0'}}
            onViewerChanged={onChartViewerChanged}
            options={spec.options}
          />,
          [[wheelGuard]],
        ) }
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
