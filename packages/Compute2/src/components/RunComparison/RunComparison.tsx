import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {COMPARISON_SERVICE, ComparisonAddition} from './comparison-service';
import {showAddToComparisonDialog} from './comparison-views';

import {IconFA, ResizeHandle, ToggleInput, Tabs, ifOverlapping} from '@datagrok-libraries/webcomponents-vue';
import {getRunTitle} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {openSaveRunDialog} from '../RunStorage/SaveRunDialog';
import {asItemMetadata, defaultRunName} from '../RunStorage/run-storage';
import {RunBrowser} from '../RunStorage/RunBrowser';

import {
  ComparisonTarget, ColumnTarget, ColumnCandidate, CandidateOverrides,
  ComparisonEntry, NameMapping, UnitMapping,
  AxisModeSelection, AxisModeConfig, AxisMode, TimeUnit, TIME_UNITS, isIndexCandidateType, isNumericType,
} from './types';
import {
  matchScalarTargets, matchColumnTargets, buildAliasGroups, buildUnitConverter,
} from './matching';
import {
  matchesFilter, unmatchedItems, isSplitCandidate, selectionToMap, computeIndexRows,
  equalTargetKeys as findEqualTargets, resolveAxisModes, targetContexts, withDefaultPicks, uniqueEntryName,
  withEntryName, withoutOutdatedOutputs,
} from './selection';
import {entryFromDataFrame, loadRunEntries} from './entry-extraction';
import {pruneOverrides, withCandidateEnabled, withoutTargetOverrides} from './target-selection';
import {columnPrimitives} from './comparison-builders';
import {showMappingsDialog} from './mappings-dialog';
import {ComparisonTab, SharedComparison, renderListFilter, sourceBadge} from './ComparisonTab';
import './RunComparison.css';

export const RunComparison = Vue.defineComponent({
  name: 'RunComparison',
  props: {
    roleOnlyFilter: {
      type: Boolean,
      default: false,
    },
  },
  setup(props) {
    const isAddingRuns = Vue.ref(false);
    const isSaving = Vue.ref(false);
    const selectedModel = Vue.shallowRef<DG.Func | undefined>(undefined);
    // models of compared snapshots that the model filter leaves out (no model role)
    const extraModels = Vue.shallowRef<DG.Func[]>([]);

    const mergeSameFuncs = Vue.ref(true);
    const hideOutdated = Vue.ref(false);

    const historySelection = Vue.shallowRef<DG.FuncCall[]>([]);
    const entries = Vue.shallowRef<ComparisonEntry[]>([]);
    // entryId -> tablePath -> index column name
    const indexSelection = Vue.ref<Record<string, Record<string, string>>>({});
    // entryId -> tablePath -> split (category) column name
    const splitSelection = Vue.ref<Record<string, Record<string, string>>>({});
    // entryId -> tablePath -> axis-mode pick; stale picks are ignored, not deleted
    const axisModeSelection = Vue.ref<AxisModeSelection>({});
    const candidateOverrides = Vue.ref<CandidateOverrides>({});

    // a tab's number names it until it charts a value
    const tabs = Vue.ref([{key: '1', number: 1}]);
    const activeTab = Vue.ref('1');
    const tabTitles = Vue.ref<Record<string, string>>({});
    let lastTabNumber = 1;

    const historyHeight = Vue.ref(440);
    const sidebarWidth = Vue.ref(360);


    const addRunsById = async (ids: string[], names = new Map<string, string>()) => {
      isAddingRuns.value = true;
      try {
        const {entries: loaded, failed} = await loadRunEntries(ids);
        for (const entry of loaded) {
          addEntry(withEntryName(entry, names.get(entry.id) ?? entry.name));
          if (entry.failedSteps)
            grok.shell.warning(`${entry.name}: could not load ${entry.failedSteps.join(', ')}`);
        }
        if (failed.length > 0)
          grok.shell.error(`Could not load ${failed.map((id) => names.get(id) ?? `run ${id}`).join(', ')}`);
      } finally {
        isAddingRuns.value = false;
      }
    };

    const add = async (
      {runIds = [], runNames, entries: added = [], model, selectModel = false}: ComparisonAddition,
    ) => {
      if (model && !extraModels.value.some((f) => f.nqName === model.nqName))
        extraModels.value = [...extraModels.value, Vue.markRaw(model)];
      if (model && selectModel)
        selectedModel.value = Vue.markRaw(model);
      const present = entries.value.filter((entry) => runIds.includes(entry.id) || runIds.includes(entry.runId ?? ''));
      if (present.length > 0)
        grok.shell.warning(`Already in this comparison: ${present.map((entry) => entry.name).join(', ')}`);
      const newRunIds = runIds.filter((id) => !hasEntry(id));
      if (newRunIds.length > 0)
        await addRunsById(newRunIds, runNames);
      added.forEach(addEntry);
    };

    const service = Vue.inject(COMPARISON_SERVICE)!;
    service.connect({add, entryNames: () => entries.value.map((entry) => entry.name)});


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

    // pre-fill pickers from the same function's picks in the other runs, then from
    // {comparisonIndex}/{comparisonSplit} output annotations; stored (user) selections always win
    const applyDefaultPicks = (entry: ComparisonEntry, others: ComparisonEntry[]) => {
      const next = withDefaultPicks(entry.nodes,
        {index: indexSelection.value, split: splitSelection.value, axis: axisModeSelection.value},
        others.map((other) => other.nodes));
      indexSelection.value = next.index;
      splitSelection.value = next.split;
      axisModeSelection.value = next.axis;
    };

    // a saved snapshot keeps its id, so its run is known by runId
    const hasEntry = (id: string) => entries.value.some((existing) => existing.id === id || existing.runId === id);

    const addEntry = (entry: ComparisonEntry) => {
      if (hasEntry(entry.id))
        return;
      const added = withEntryName(entry, uniqueEntryName(entry.name, entries.value.map((existing) => existing.name)));
      applyDefaultPicks(added, entries.value);
      entries.value = [...entries.value, Vue.markRaw(added)];
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

    const addSelectedRuns = () => {
      if (isAddingRuns.value || historySelection.value.length === 0)
        return;
      const calls = historySelection.value.filter((call) => !hasEntry(call.id));
      if (calls.length > 0)
        showAddToComparisonDialog(calls.map((call) => ({runId: call.id, name: getRunTitle(call)})), {into: service});
    };

    const saveSnapshot = async (entry: ComparisonEntry) => {
      const {metadata, model} = entry.snapshot!;
      if (!model)
        return;
      const result = await openSaveRunDialog({func: model, name: defaultRunName(model, undefined, metadata.title),
        description: metadata.description, tags: metadata.tags});
      if (!result)
        return;
      isSaving.value = true;
      try {
        const run = await entry.snapshot!.save!(asItemMetadata(result));
        const others = entries.value.filter((item) => item.id !== entry.id).map((item) => item.name);
        const named = withEntryName(entry, uniqueEntryName(getRunTitle(run), others));
        const saved: ComparisonEntry = {...named, runId: run.id, snapshot: undefined};
        entries.value = entries.value.map((item) => item.id === entry.id ? Vue.markRaw(saved) : item);
      } catch (e: any) {
        grok.shell.error(e);
      } finally {
        isSaving.value = false;
      }
    };

    const removeWithConfirm = (entry: ComparisonEntry) => {
      if (!entry.snapshot) {
        removeEntry(entry.id);
        return;
      }
      ui.dialog('Remove confirmation')
        .add(ui.divText(`'${entry.name}' is not stored in history, so it cannot be added back. Remove it?`))
        .onOK(() => removeEntry(entry.id))
        .show({center: true, modal: true});
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

    const hasOutdated = Vue.computed(() => entries.value.some((entry) => entry.outdatedSteps));
    const entryNodes = Vue.computed(() => entries.value.map((entry) =>
      hideOutdated.value ? withoutOutdatedOutputs(entry.nodes) : entry.nodes));

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
      const nodes = entryNodes.value;
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
      candidateOverrides.value = withCandidateEnabled(candidateOverrides.value, target, candidate, enabled);
    };

    const resetTargetOverrides = (targetKey: string) => {
      candidateOverrides.value = withoutTargetOverrides(candidateOverrides.value, targetKey);
    };

    // removed entries are not consumed by matching, so the override rewrite converges in one cycle
    Vue.watch(targets, (list) => {
      candidateOverrides.value = pruneOverrides(candidateOverrides.value, list);
    });

    const indexFilter = Vue.ref('');

    const columnValues = (entryId: string, tablePath: string, columnName: string) => {
      const col = entriesById.value.get(entryId)?.dataFrames.get(tablePath)?.col(columnName);
      return col ? columnPrimitives(col) : undefined;
    };

    const equalKeys = Vue.computed(() => findEqualTargets(targets.value, columnValues, indexColumnType));

    const editMappings = () => showMappingsDialog(nameMappings.value, unitMappings.value, unmatchedItems(
      entryNodes.value, targets.value, indexColumnsMap.value,
    ), (names, units) => {
      nameMappings.value = names;
      unitMappings.value = units;
    });

    // non-default axis modes only; tables whose index moved off a time type drop out
    const resolvedAxisModes = Vue.computed(() => resolveAxisModes(
      entryNodes.value,
      indexColumnsMap.value,
      axisModeSelection.value,
    ));

    const shared: SharedComparison = {
      entries, entriesById, entryNodes, targets, contexts, equalKeys, indexColumnsMap, resolvedAxisModes,
      aliasGroups, unitConverter, candidateOverrides, mappingCount, hasOutdated, hideOutdated,
      indexColumnType, setCandidateEnabled, resetTargetOverrides, editMappings,
    };

    const addTab = () => {
      const number = ++lastTabNumber;
      tabs.value = [...tabs.value, {key: String(number), number}];
      activeTab.value = String(number);
    };

    // closing the shown tab shows its right neighbor, or the left one for the last tab
    const closeTab = (key: string) => {
      const index = tabs.value.findIndex((tab) => tab.key === key);
      const rest = tabs.value.filter((tab) => tab.key !== key);
      if (activeTab.value === key)
        activeTab.value = rest[Math.min(index, rest.length - 1)].key;
      tabs.value = rest;
      const {[key]: _removed, ...titles} = tabTitles.value;
      tabTitles.value = titles;
    };

    const tabItems = Vue.computed(() => tabs.value.map((tab) => ({
      key: tab.key,
      title: tabTitles.value[tab.key] || `Comparison ${tab.number}`,
      closable: tabs.value.length > 1,
    })));

    // tables that could participate in column comparison and need an index choice;
    // with merging on, same-function outputs (by nqName) collapse into one row
    const indexRows = Vue.computed(() => computeIndexRows(
      entryNodes.value,
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
              { entry.snapshot &&
                <span class='c2-comparison-badge c2-badge-unsaved' title='Current state, not stored in history'>
                  unsaved
                </span> }
              <span style={{overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap', flex: '1'}}
                title={entry.modelName ? `${entry.modelName}: ${entry.name}` : entry.name}>
                {entry.name}
              </span>
              { entry.outdatedSteps &&
                <IconFA name='clock' tooltip={`Outdated outputs: ${entry.outdatedSteps.join(', ')}`}/> }
              { entry.snapshot?.save &&
                <IconFA name='save' tooltip='Save to history' onClick={() => saveSnapshot(entry)}/> }
              <IconFA name='times' tooltip='Remove from comparison' onClick={() => removeWithConfirm(entry)}/>
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

    return () => Vue.withDirectives(
      <div style={{display: 'flex', height: '100%', width: '100%', overflow: 'hidden'}}>
        <div style={{
          width: `${sidebarWidth.value}px`, flexShrink: '0', display: 'flex', flexDirection: 'column',
          gap: '8px', padding: '8px', overflow: 'auto',
        }}>
          <div style={{fontWeight: 'bold'}}>Model</div>
          <div style={{
            height: `${historyHeight.value}px`, flexShrink: '0',
            border: '1px solid var(--grey-1)', borderRadius: '3px',
          }}>
            <RunBrowser
              model={selectedModel.value}
              onUpdate:model={(func?: DG.Func) => selectedModel.value = func && Vue.markRaw(func)}
              selectModel={true}
              allowAllModels={true}
              extraModels={extraModels.value}
              roleOnlyFilter={props.roleOnlyFilter}
              allowCompare={false}
              allowOpen={false}
              showDetails={false}
              onSelectionChanged={(calls: DG.FuncCall[]) => historySelection.value = calls.map((c) => Vue.markRaw(c))}
            />
          </div>
          <ResizeHandle
            axis='y'
            size={historyHeight.value}
            min={150}
            max={800}
            onUpdate:size={(size) => historyHeight.value = size}
          />
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
          { entries.value.length < 2 &&
            <div style={{color: 'var(--grey-4)'}}>Add at least two runs (or tables) to compare</div> }
          { /* hidden, not unmounted, so the tabs keep their state while runs are swapped */ }
          <div style={{
            display: entries.value.length < 2 ? 'none' : 'flex',
            flexDirection: 'column', flex: '1', minHeight: '0px', gap: '12px',
          }}>
            { renderIndexPickers() }
            <Tabs
              items={tabItems.value}
              activeKey={activeTab.value}
              addable={true}
              style={{flexShrink: '1', minHeight: '300px'}}
              onUpdate:activeKey={(key: string) => activeTab.value = key}
              onClose={closeTab}
              onAdd={addTab}
            >
              {{pane: ({item}: {item: {key: string}}) =>
                <ComparisonTab
                  shared={shared}
                  onTitle={(title: string) => tabTitles.value = {...tabTitles.value, [item.key]: title}}
                />}}
            </Tabs>
          </div>
        </div>
      </div>,
      [[ifOverlapping, isAddingRuns.value || isSaving.value, 'Loading...']],
    );
  },
});
