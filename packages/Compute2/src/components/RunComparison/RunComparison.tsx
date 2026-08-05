import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {onClickOutside} from '@vueuse/core';

import {
  IconFA, Viewer, ResizeHandle, ToggleInput, ifOverlapping, wheelGuard,
} from '@datagrok-libraries/webcomponents-vue';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {getModelFilter} from '@datagrok-libraries/compute-utils/model-catalog/src/model-handler';

import {History} from '../History/History';
import {
  ComparisonTarget, ScalarTarget, ColumnTarget, ColumnCandidate, CandidateOverrides, MatchConfidence,
  ComparisonEntry, EntrySourceKind, RUN_COLUMN, candidateId,
  AxisModeSelection, AxisModeConfig, AxisMode, TimeUnit, TIME_UNITS, isIndexCandidateType,
} from './types';
import {matchScalarTargets, matchColumnTargets} from './matching';
import {
  getEntryStatuses, matchesFilter, compatibleTargetsFor, multiValueOverlap,
  isSplitCandidate, selectionToMap, computeIndexRows, isTargetEqualAcrossRuns, resolveAxisModes,
} from './selection';
import {entryFromFuncCall, entryFromDataFrame} from './entry-extraction';
import {
  buildScalarComparison, buildMultiScalarComparison, buildColumnComparison, buildMultiColumnComparison,
} from './comparison-builders';
import './RunComparison.css';

// chip styling lives in RunComparison.css (c2-comparison-badge + these modifiers)
const sourceBadge = (kind: EntrySourceKind) =>
  <span class={`c2-comparison-badge c2-badge-${kind}`}>{kind}</span>;

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
    const selectedTargetKey = Vue.ref('');
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
    tableInput.onChanged.subscribe(() => {
      const df = tableInput.value;
      if (!df)
        return;
      addEntry(entryFromDataFrame(df));
      // resetting synchronously re-enters the Dart stream controller mid-dispatch
      setTimeout(() => tableInput.value = null);
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

    const indexColumnType = (entryId: string, tablePath: string, columnName: string) => {
      const entry = entries.value.find((item) => item.id === entryId);
      const table = entry?.nodes.tables.find((item) => item.path === tablePath);
      return table?.columns.find((col) => col.name === columnName)?.type;
    };

    const defaultIndexOf = (entryId: string, tablePath: string) => {
      const entry = entries.value.find((item) => item.id === entryId);
      return entry?.nodes.tables.find((item) => item.path === tablePath)?.defaultIndexColumn;
    };

    const indexColumnsMap = Vue.computed(() => selectionToMap(indexSelection.value,
      (entryId, tablePath, columnName) =>
        isIndexCandidateType(indexColumnType(entryId, tablePath, columnName)) ||
        columnName === defaultIndexOf(entryId, tablePath)));

    const splitColumnsMap = Vue.computed(() => selectionToMap(splitSelection.value,
      (entryId, tablePath, columnName) => isSplitCandidate(
        {name: columnName, type: indexColumnType(entryId, tablePath, columnName) ?? ''},
        indexColumnsMap.value.get(entryId)?.get(tablePath) ?? '')));

    const targets = Vue.computed<ComparisonTarget[]>(() => {
      const nodes = entries.value.map((entry) => entry.nodes);
      return [
        ...matchScalarTargets(nodes),
        ...matchColumnTargets(nodes, indexColumnsMap.value, splitColumnsMap.value, candidateOverrides.value),
      ].sort((a, b) => a.kind !== b.kind ?
        (a.kind === 'scalar' ? -1 : 1) :
        a.displayName.localeCompare(b.displayName));
    });

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

    // drop overrides that no longer resolve (entry removed, index changed, cluster reshaped);
    // removed entries are not consumed by matching, so the rewrite converges in one cycle
    Vue.watch(targets, (list) => {
      const validIds = new Map(list
        .filter((target): target is ColumnTarget => target.kind === 'column')
        .map((target) => [target.key,
          new Set(target.candidates.map((candidate) => candidateId(candidate.binding)))]));
      let changed = false;
      const next: CandidateOverrides = {};
      for (const [key, ids] of Object.entries(candidateOverrides.value)) {
        const valid = validIds.get(key);
        const kept = valid ? Object.entries(ids).filter(([id]) => valid.has(id)) : [];
        if (kept.length !== Object.keys(ids).length)
          changed = true;
        if (kept.length > 0)
          next[key] = Object.fromEntries(kept);
      }
      if (changed)
        candidateOverrides.value = next;
      // drop multi-value selections whose target vanished (structural change);
      // merely-conflicting targets keep their keys so reverting an edit restores them
      const keys = new Set(list.map((target) => target.key));
      const keptKeys = multiKeys.value.filter((key) => keys.has(key));
      if (keptKeys.length !== multiKeys.value.length) {
        multiKeys.value = keptKeys;
        if (keptKeys.length === 0)
          multiMode.value = false;
      }
    });

    const indexFilter = Vue.ref('');
    const targetFilter = Vue.ref('');

    // dayjs and other value objects collapse to primitives so cross-run compares work
    const columnValues = (entryId: string, tablePath: string, columnName: string) => {
      const entry = entries.value.find((item) => item.id === entryId);
      const col = entry?.dataFrames.get(tablePath)?.col(columnName);
      if (!col)
        return undefined;
      const values = new Array(col.length);
      for (let i = 0; i < col.length; i++) {
        const value = col.isNone(i) ? null : col.get(i);
        values[i] = value instanceof Object ? value.valueOf() : value;
      }
      return values;
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
        .filter((target) => isTargetEqualAcrossRuns(target, cached))
        .map((target) => target.key));
    });

    // in multi mode: suggestions plus everything already selected — a selected target
    // that became conflicting must stay visible so it can be unchecked
    const filteredTargets = Vue.computed(() => {
      const listed: ComparisonTarget[] = multiMode.value ?
        targets.value.filter((target) =>
          multiKeys.value.includes(target.key) ||
            compatibleTargets.value.some((item) => item.key === target.key)) :
        targets.value;
      return listed.filter((target) =>
        (!equalTargetKeys.value.has(target.key) ||
          (multiMode.value && multiKeys.value.includes(target.key))) &&
        matchesFilter(targetFilter.value, target.displayName));
    });

    const selectedTarget = Vue.computed(() =>
      targets.value.find((target) => target.key === selectedTargetKey.value) ?? null);

    const multiMode = Vue.ref(false);
    const multiKeys = Vue.ref<string[]>([]);

    const scalarChartType = Vue.ref<'radar' | 'pcplot'>('radar');
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

    // snapshot export: clone of the chart data plus the chart with its current options
    const snapshotView = (name: string, addToWorkspace: boolean) => {
      const df = comparison.value!.chartDf.clone();
      df.name = name;
      const view = addToWorkspace ? grok.shell.addTableView(df) : DG.TableView.create(df, false);
      const options = chartViewer.value?.getOptions();
      if (options)
        view.addViewer(options.type, options.look);
      return view;
    };

    // the platform Save-project dialog handles upload, layout linking, and sharing.
    // Project.showSaveDialog is newer than the last released js-api (added 2026-07-16,
    // after datagrok-api 1.27.7), so it is reached via a cast and feature-detected at
    // runtime; older platforms may have the Dart binding but not the wrapper
    const jsApiShowSaveDialog = (DG.Project as any).showSaveDialog;
    const canSaveProject = typeof jsApiShowSaveDialog === 'function' ||
      typeof (window as any).grok_Project_OpenSaveDialog === 'function';

    const showProjectSaveDialog = (view: DG.TableView, name: string): Promise<DG.Project | null> => {
      if (typeof jsApiShowSaveDialog === 'function')
        return jsApiShowSaveDialog.call(DG.Project, {tables: [view.dataFrame], views: [view], name});
      return (window as any).grok_Project_OpenSaveDialog(
        [view.dataFrame.dart], [view.dart], [], name, '', '');
    };

    const saveAndShare = async (name: string) => {
      const view = snapshotView(name, false);
      try {
        await showProjectSaveDialog(view, name);
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        view.detach();
      }
    };

    const exportComparison = () => {
      const currentComparison = comparison.value;
      if (!currentComparison || currentComparison.chartDf.rowCount === 0)
        return;
      const nameInput = ui.input.string('Name', {value: currentComparison.chartDf.name || 'Comparison'});
      const dlg = ui.dialog('Export comparison').add(nameInput.root);
      const openLabel = 'OPEN IN WORKSPACE';
      if (canSaveProject) {
        dlg.addButton(openLabel, () => {
          dlg.close();
          snapshotView(nameInput.value || 'Comparison', true);
        });
        dlg.onOK(() => saveAndShare(nameInput.value || 'Comparison'));
        dlg.show({center: true});
        dlg.getButton('OK').innerText = 'SAVE & SHARE';
      } else {
        dlg.onOK(() => snapshotView(nameInput.value || 'Comparison', true));
        dlg.show({center: true});
        dlg.getButton('OK').innerText = openLabel;
      }
    };

    // DG inputs fire onChanged on programmatic sets too, so rendering a new mode value
    // echoes back through the toggle — a no-op update must not clobber multiKeys
    const setMultiMode = (val: boolean) => {
      if (val === multiMode.value)
        return;
      if (val)
        multiKeys.value = selectedTargetKey.value ? [selectedTargetKey.value] : [];
      else if (multiKeys.value.length > 0)
        selectedTargetKey.value = multiKeys.value[0];
      multiMode.value = val;
    };

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
      resolvedAxisModes.value,
    ));

    const comparison = Vue.computed(() => {
      const target = selectedTarget.value;
      if (!target || entries.value.length < 2)
        return null;
      if (target.kind === 'scalar') {
        if (multiMode.value) {
          const selected = targets.value.filter((item): item is ScalarTarget =>
            item.kind === 'scalar' && multiKeys.value.includes(item.key));
          const anchorIndex = selected.findIndex((item) => item.key === target.key);
          if (anchorIndex > 0)
            selected.unshift(...selected.splice(anchorIndex, 1));
          if (selected.length === 0)
            return null;
          if (selected.length > 1) {
            const result = buildMultiScalarComparison(selected, entries.value);
            return result ? Vue.markRaw({kind: 'multi-scalar' as const, target, ...result}) : null;
          }
          const result = buildScalarComparison(selected[0], entries.value);
          return Vue.markRaw({kind: 'scalar' as const, target: selected[0], ...result});
        }
        const result = buildScalarComparison(target, entries.value);
        return Vue.markRaw({kind: 'scalar' as const, target, ...result});
      }
      if (multiMode.value) {
        const selected = targets.value.filter((item): item is ColumnTarget =>
          item.kind === 'column' && multiKeys.value.includes(item.key));
        const anchorIndex = selected.findIndex((item) => item.key === target.key);
        if (anchorIndex > 0)
          selected.unshift(...selected.splice(anchorIndex, 1));
        if (selected.length === 0)
          return null;
        if (selected.length > 1) {
          const result = buildMultiColumnComparison(selected, entries.value, resolvedAxisModes.value);
          return result ? Vue.markRaw({kind: 'column' as const, target, ...result}) : null;
        }
        const result = buildColumnComparison(selected[0], entries.value, resolvedAxisModes.value);
        return result ? Vue.markRaw({kind: 'column' as const, target: selected[0], ...result}) : null;
      }
      const result = buildColumnComparison(target, entries.value, resolvedAxisModes.value);
      return result ? Vue.markRaw({kind: 'column' as const, target, ...result}) : null;
    });

    // a 2-axis radar degenerates to a line; force the PC plot (a slopegraph) until a 3rd value is added
    const radarEligible = Vue.computed(() => {
      const current = comparison.value;
      return radarAvailable && current?.kind === 'multi-scalar' && current.valueColumnNames.length >= 3;
    });
    const effectiveScalarChartType = Vue.computed(() =>
      radarEligible.value ? scalarChartType.value : 'pcplot');

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
        const entry = entries.value.find((item) => item.id === candidate.binding.entryId);
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
            const isExpanded = target.kind === 'column' && !!expandedTargetKeys.value[target.key];
            const anchor = selectedTarget.value;
            const overlap = multiMode.value && anchor && target.kind === anchor.kind &&
              target.key !== anchor.key ?
              multiValueOverlap(anchor, target) : null;
            const gapRuns = overlap ? [...overlap.missing, ...overlap.conflicting]
              .map((id) => entries.value.find((entry) => entry.id === id)?.name ?? id) : [];
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
              <IconFA
                name={isExpanded ? 'chevron-down' : 'chevron-right'}
                tooltip='Show matched items'
                onClick={(e: Event) => {e.stopPropagation(); toggleExpanded(target.key);}}
              /> :
              <span></span> }
            <IconFA name={target.kind === 'scalar' ? 'hashtag' : 'table'} tooltip={target.kind}/>
            <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap'}}
              title={target.displayName}>
              {target.displayName}
            </span>
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
                title={anchor!.kind === 'scalar' ?
                  `Missing value on: ${gapRuns.join(', ')}` :
                  `Shown as gaps (no shared rows with ${anchor!.displayName}): ${gapRuns.join(', ')}`}
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
          { (multiMode.value || compatibleTargets.value.length > 1) &&
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
                cursor: 'pointer', color: 'var(--blue-1)',
              }}
              onClick={exportComparison}
            >
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
          const entry = entries.value.find((item) => item.id === status.entryId);
          if (!entry)
            return null;
          if (!status.matched)
            return {entry, text: status.reason, warning: false};
          if (status.warning)
            return {entry, text: status.warning, warning: true};
          return null;
        })
        .filter((item) => item != null);
      if (chips.length === 0)
        return null;
      return <div style={{display: 'flex', gap: '6px', flexWrap: 'wrap', padding: '4px 0px'}}>
        { chips.map(({entry, text, warning}) =>
          <span
            key={entry.id}
            class={`c2-comparison-status ${warning ? 'c2-status-warning' : 'c2-status-error'}`}
            title={warning ?
              `The value charts as a plain series until every participating table has the ` +
              `${text === 'independent points not set' ? 'independent points' : 'relative timeseries'} mode set` :
              entry.name}>
            {entry.name}: {text}
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
      // multi-value needs a line-chartable index (suggestions are gated on it), but an index
      // switched to a string column after selection would land the frame on a line chart
      if (currentComparison.kind === 'column' && currentComparison.isKeyIndex &&
        'valueColumnNames' in currentComparison) {
        return <div style={{color: 'var(--grey-4)', padding: '8px 0px'}}>
          Multiple values cannot chart over a categorical index — pick a numeric or datetime
          index, or select a single value
        </div>;
      }
      const chartMinHeight = currentComparison.kind === 'multi-scalar' ? 300 :
        currentComparison.kind === 'column' && currentComparison.isScatter ? 300 :
          'valueColumnNames' in currentComparison ?
            250 * Math.max(1, (currentComparison.valueColumnNames as string[]).length) : 250;
      const effectiveChartHeight = Math.max(chartMinHeight, chartHeight.value);
      const chartStyle = {width: '100%', height: `${effectiveChartHeight}px`, flexShrink: '0'};
      let chart;
      if (currentComparison.kind === 'multi-scalar') {
        chart = effectiveScalarChartType.value === 'radar' ?
          <Viewer
            type='Radar'
            dataFrame={currentComparison.chartDf}
            style={chartStyle}
            onViewerChanged={onChartViewerChanged}
            options={{
              valuesColumnNames: currentComparison.valueColumnNames,
              colorColumnName: RUN_COLUMN,
              showCurrentRow: false,
            }}
          /> :
          <Viewer
            type={DG.VIEWER.PC_PLOT}
            dataFrame={currentComparison.chartDf}
            style={chartStyle}
            onViewerChanged={onChartViewerChanged}
            options={{
              columnNames: currentComparison.valueColumnNames,
              colorColumnName: RUN_COLUMN,
            }}
          />;
      } else if (currentComparison.kind === 'scalar') {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={{
            valueColumnName: currentComparison.valueColumnName,
            valueAggrType: 'avg',
            splitColumnName: RUN_COLUMN,
          }}
        />;
      } else if (currentComparison.isKeyIndex) {
        chart = <Viewer
          type={DG.VIEWER.BAR_CHART}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={{
            valueColumnName: currentComparison.valueColumnName,
            valueAggrType: 'avg',
            splitColumnName: currentComparison.indexColumnName,
            stackColumnName: RUN_COLUMN,
          }}
        />;
      } else if (currentComparison.isScatter) {
        // color: value name (multi-value melt) > split > Run; run identity moves to
        // marker shapes whenever color is taken
        const melted = currentComparison.melted;
        const splitColumnName = currentComparison.splitColumnName;
        chart = <Viewer
          type={DG.VIEWER.SCATTER_PLOT}
          dataFrame={currentComparison.chartDf}
          style={chartStyle}
          onViewerChanged={onChartViewerChanged}
          options={ melted ? {
            xColumnName: currentComparison.indexColumnName,
            yColumnName: melted.valueColumnName,
            colorColumnName: melted.seriesColumnName,
            markersColumnName: RUN_COLUMN,
          } : {
            xColumnName: currentComparison.indexColumnName,
            yColumnName: currentComparison.valueColumnName,
            colorColumnName: splitColumnName ?? RUN_COLUMN,
            ...splitColumnName ? {markersColumnName: RUN_COLUMN} : {},
          }}
        />;
      } else {
        const yColumnNames = 'valueColumnNames' in currentComparison ?
          currentComparison.valueColumnNames as string[] : [currentComparison.valueColumnName];
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
        { currentComparison.kind === 'multi-scalar' && radarEligible.value &&
          <div style={{display: 'flex', gap: '4px', padding: '2px 0px', flexShrink: '0'}}>
            { ([['radar', 'Radar'], ['pcplot', 'PC Plot']] as const).map(([value, label]) =>
              <span
                key={value}
                onClick={() => scalarChartType.value = value}
                style={{
                  cursor: 'pointer', fontSize: '11px', borderRadius: '3px', padding: '1px 8px',
                  border: '1px solid var(--grey-2)',
                  background: effectiveScalarChartType.value === value ?
                    'var(--blue-1)' : 'transparent',
                  color: effectiveScalarChartType.value === value ? 'var(--white)' : 'var(--grey-5)',
                }}
              >{label}</span>,
            )}
          </div> }
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
