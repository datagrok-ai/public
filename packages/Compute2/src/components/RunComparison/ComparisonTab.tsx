import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';

import {IconFA, Viewer, ResizeHandle, ToggleInput, wheelGuard} from '@datagrok-libraries/webcomponents-vue';
import {
  AxisConfigMap, CandidateOverrides, ColumnCandidate, ColumnTarget, ComparisonEntry, ComparisonEntryNodes,
  ComparisonTarget, EntrySourceKind, MatchConfidence, candidateId,
} from './types';
import {UnitConverter} from './matching';
import {
  compatibleTargetsFor, getEntryStatuses, groupCandidatesByEntry, listedTargets, selectionGaps, statusChips,
} from './selection';
import {
  TargetSelection, EMPTY_TARGET_SELECTION, isCharted, clickTarget, setMultiMode, reconcileTargets, chartedTargets,
} from './target-selection';
import {canSaveProject, captureTableLayout, showProjectSaveDialog} from '../../project-export';
import {buildComparison} from './comparison-builders';
import {ScalarChartType, chartSpec} from './chart-spec';

// chip styling lives in RunComparison.css (c2-comparison-badge + these modifiers)
export const sourceBadge = (kind: EntrySourceKind) =>
  <span class={`c2-comparison-badge c2-badge-${kind}`}>{kind}</span>;

const confidenceBadge = (confidence: MatchConfidence) =>
  <span class={`c2-comparison-badge c2-confidence-${confidence}`}>{confidence}</span>;

export const renderListFilter = (value: Vue.Ref<string>, placeholder: string) => (
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

type ColumnTypeGetter = (entryId: string, tablePath: string, columnName: string) => string | undefined;

/** What every tab of a view shares: the compared runs, their matched values, picks and mappings. */
export interface SharedComparison {
  entries: Readonly<Vue.Ref<ComparisonEntry[]>>;
  entriesById: Readonly<Vue.Ref<Map<string, ComparisonEntry>>>;
  entryNodes: Readonly<Vue.Ref<ComparisonEntryNodes[]>>;
  targets: Readonly<Vue.Ref<ComparisonTarget[]>>;
  contexts: Readonly<Vue.Ref<Map<string, string>>>;
  // targets equal across runs, for the tabs that hide them
  equalKeys: Readonly<Vue.Ref<Set<string>>>;
  indexColumnsMap: Readonly<Vue.Ref<Map<string, Map<string, string>>>>;
  resolvedAxisModes: Readonly<Vue.Ref<AxisConfigMap>>;
  aliasGroups: Readonly<Vue.Ref<Map<string, number>>>;
  unitConverter: Readonly<Vue.Ref<UnitConverter>>;
  candidateOverrides: Readonly<Vue.Ref<CandidateOverrides>>;
  mappingCount: Readonly<Vue.Ref<number>>;
  hasOutdated: Readonly<Vue.Ref<boolean>>;
  hideOutdated: Vue.Ref<boolean>;
  indexColumnType: ColumnTypeGetter;
  setCandidateEnabled: (target: ColumnTarget, candidate: ColumnCandidate, enabled: boolean) => void;
  resetTargetOverrides: (targetKey: string) => void;
  editMappings: () => void;
}

/** One comparison of the shared runs: the Compare list with its own selection, and the chart. */
export const ComparisonTab = Vue.defineComponent({
  name: 'ComparisonTab',
  props: {
    shared: {
      type: Object as Vue.PropType<SharedComparison>,
      required: true,
    },
  },
  emits: {
    // the charted value names the tab; empty when nothing is charted
    title: (_title: string) => true,
  },
  setup(props, {emit}) {
    const {
      entries, entriesById, entryNodes, targets, contexts, equalKeys, indexColumnsMap, resolvedAxisModes,
      aliasGroups, unitConverter, candidateOverrides, mappingCount, hasOutdated, hideOutdated,
      indexColumnType, setCandidateEnabled, resetTargetOverrides, editMappings,
    } = props.shared;

    const hideEqual = Vue.ref(true);
    const targetSelection = Vue.shallowRef<TargetSelection>(EMPTY_TARGET_SELECTION);
    const selectedTargetKey = Vue.computed(() => targetSelection.value.selectedKey);
    const multiMode = Vue.computed(() => targetSelection.value.multiMode);
    const expandedTargetKeys = Vue.ref<Record<string, boolean>>({});
    const targetFilter = Vue.ref('');
    // user override for the results chart height; the per-type minimum still applies
    const chartHeight = Vue.ref(0);

    const toggleExpanded = (targetKey: string) => {
      expandedTargetKeys.value = {
        ...expandedTargetKeys.value,
        [targetKey]: !expandedTargetKeys.value[targetKey],
      };
    };

    Vue.watch(targets, (list) => {
      targetSelection.value = reconcileTargets(targetSelection.value, new Set(list.map((target) => target.key)));
    });

    const equalTargetKeys = Vue.computed(() => hideEqual.value ? equalKeys.value : new Set<string>());

    const filteredTargets = Vue.computed(() => listedTargets(targets.value, {
      multiMode: multiMode.value,
      isCharted: (key) => isCharted(targetSelection.value, key),
      compatible: compatibleTargets.value,
      equalKeys: equalTargetKeys.value,
      filter: targetFilter.value,
      contexts: contexts.value,
    }));

    const selectedTarget = Vue.computed(() =>
      targets.value.find((target) => target.key === selectedTargetKey.value) ?? null);

    const scalarChartType = Vue.ref<ScalarChartType>('radar');
    const sharedScale = Vue.ref(true);
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

    // DG inputs fire onChanged on programmatic sets too, so rendering a new mode value
    // echoes back through the toggle — the no-op transition keeps multiKeys intact
    const onMultiModeChanged = (on: boolean) => {
      targetSelection.value = setMultiMode(targetSelection.value, on);
    };

    const compatibleKeys = (anchorKey: string) => compatibleTargetsFor(
      targets.value.find((target) => target.key === anchorKey) ?? null, targets.value, indexColumnType,
    ).map((target) => target.key);

    const entryStatuses = Vue.computed(() => getEntryStatuses(
      entryNodes.value,
      selectedTarget.value,
      indexColumnsMap.value,
      {axisModes: resolvedAxisModes.value, aliases: aliasGroups.value, convert: unitConverter.value},
    ));

    const charted = Vue.computed(() => chartedTargets(targetSelection.value, targets.value));

    const gaps = Vue.computed(() => multiMode.value ? selectionGaps(charted.value) : new Map<string, string[]>());

    const comparison = Vue.computed(() => {
      const result = buildComparison(charted.value, entries.value, resolvedAxisModes.value);
      return result ? Vue.markRaw(result) : null;
    });

    const currentChartSpec = Vue.computed(() =>
      chartSpec(comparison.value, selectedTarget.value != null, scalarChartType.value, radarAvailable,
        sharedScale.value));

    const title = Vue.computed(() => {
      const [first, ...rest] = charted.value;
      return first ? `${first.displayName}${rest.length > 0 ? ` +${rest.length}` : ''}` : '';
    });
    Vue.watch(title, (value) => emit('title', value), {immediate: true});

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
        { candidate.unitsAssumed &&
          <IconFA name='info-circle' tooltip={`Raw data without units is taken to be in ${target.units}`}/> }
        { candidate.auto &&
          <span style={{fontSize: '11px', color: 'var(--grey-4)', flexShrink: '0'}}>(auto)</span> }
      </label>;
    };

    const renderCandidatePanel = (target: ColumnTarget) => {
      const groups = groupCandidatesByEntry(target.candidates).flatMap(({entryId, candidates}) => {
        const entry = entriesById.value.get(entryId);
        return entry ? [{entry, candidates}] : [];
      });
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
            { hasOutdated.value &&
              <ToggleInput
                caption='Hide outdated outputs'
                tooltip='Leave out outputs of steps whose inputs changed after the last run, or that never ran'
                value={hideOutdated.value}
                onUpdate:value={(val) => hideOutdated.value = val}
              /> }
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
            gridTemplateColumns:
              'max-content max-content fit-content(400px) max-content max-content max-content max-content 1fr',
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
              onMousedown={(e: MouseEvent) => {if (e.shiftKey) e.preventDefault();}}
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
                target.kind === 'column' && target.rawWithoutUnits ?
                  <IconFA name='info-circle' tooltip={`Raw data without units is taken to be in ${target.units}`}/> :
                  <span></span> }
              { confidenceBadge(target.confidence) }
              <span style={{fontSize: '11px', color: 'var(--grey-4)', flexShrink: '0'}}>
                {target.coverage}/{target.total}
              </span>
              { gapRuns.length > 0 ?
                <span
                  class='c2-comparison-badge c2-badge-partial'
                  title={`No value on: ${gapRuns.join(', ')}`}
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
      const chips = statusChips(entryStatuses.value).flatMap((chip) => {
        const entry = entriesById.value.get(chip.entryId);
        return entry ? [{...chip, entry}] : [];
      });
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
        { spec.sharedScaleSwitch &&
          <div style={{padding: '2px 0px', flexShrink: '0'}}>
            <ToggleInput
              caption='Shared y axis'
              value={sharedScale.value}
              onUpdate:value={(val) => sharedScale.value = val}
            />
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

    return () => (
      <div style={{
        display: 'flex', flexDirection: 'column', flex: '1', minHeight: '0px', gap: '12px', paddingTop: '8px',
      }}>
        { renderTargets() }
        { renderResultsHeader() }
        { renderStatuses() }
        { renderComparison() }
      </div>
    );
  },
});
