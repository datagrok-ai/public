import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {PipelineState, isFuncCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Button} from '@datagrok-libraries/webcomponents-vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {Logger} from '../Logger/Logger';
import {PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import VueJsonPretty from 'vue-json-pretty';
import 'vue-json-pretty/lib/styles.css';
import {LinksData} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {MatchedNodePaths} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/link-matching';
import {formatNodePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FuncCallStateInfo, ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {LinkIOParsed, LinkSelectorSegment, LinkTagSegment} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/LinkSpec';
import {BehaviorSubject} from 'rxjs';
import {FilterDropdown, FilterOption} from './FilterDropdown';

interface StepStates {
  calls: Record<string, FuncCallStateInfo | undefined>;
  validations: Record<string, Record<string, ValidationResult> | undefined>;
  consistency: Record<string, Record<string, ConsistencyInfo> | undefined>;
  meta: Record<string, Record<string, BehaviorSubject<any>> | undefined>;
  descriptions: Record<string, Record<string, string | string[]> | undefined>;
}

// ---- Serialization helpers ----

function summarizeDataFrame(df: DG.DataFrame): string {
  return `#DataFrame(${df.rowCount} rows, ${df.columns.length} cols)`;
}

function summarizeValue(v: any): any {
  return v instanceof DG.DataFrame ? summarizeDataFrame(v) : v;
}

function summarizeFuncCallParams(params: DG.FuncCallParam[]): Record<string, any> {
  return Object.fromEntries(params.map((p) => [p.name, summarizeValue(p.value)]));
}

/** Check if an array looks like LinkIOParsed[] (has name + segments on first element). */
function isLinkIOParsedArray(v: any): v is LinkIOParsed[] {
  return Array.isArray(v) && v.length > 0
    && typeof v[0]?.name === 'string'
    && Array.isArray(v[0]?.segments);
}

/** Reconstruct a spec string from a parsed LinkIOParsed object.
 *  e.g. {name: "in1", segments: [{selector: "first", ids: ["stepAdd"]}, ...]} → "in1:stepAdd/res" */
function linkIOToString(io: LinkIOParsed): string {
  const parts: string[] = [io.name];
  if (io.flags?.length)
    parts[0] += `(${io.flags.join(',')})`;
  for (const seg of io.segments) {
    if (seg.type === 'tag') {
      const t = seg as LinkTagSegment;
      const ref = t.ref ? `@${t.ref},` : '';
      parts.push(`#${t.selector}(${ref}${t.tags.join('&')})`);
    } else {
      const s = seg as LinkSelectorSegment;
      if (s.selector === 'first' && !s.ref && s.stopIds.length === 0)
        parts.push(s.ids.join('|'));
      else {
        const ref = s.ref ? `@${s.ref},` : '';
        const stop = s.stopIds.length ? `,${s.stopIds.join('|')}` : '';
        parts.push(`${s.selector}(${ref}${s.ids.join('|')}${stop})`);
      }
    }
  }
  return parts[0] + (parts.length > 1 ? ':' + parts.slice(1).join('/') : '');
}

/** Convert LinkIOParsed[] to readable spec strings. */
function linkIOArrayToStrings(ios: LinkIOParsed[] | undefined): string | string[] | undefined {
  if (!ios || ios.length === 0) return undefined;
  const strs = ios.map(linkIOToString);
  return strs.length === 1 ? strs[0] : strs;
}

/** Convert a MatchedNodePaths record to readable path strings. */
function matchedPathsToStrings(paths: Record<string, MatchedNodePaths>): Record<string, string[]> {
  const result: Record<string, string[]> = {};
  for (const [name, matched] of Object.entries(paths))
    result[name] = matched.map((m) => formatNodePath(m.path, m.ioName));
  return result;
}

/** JSON replacer that handles all Datagrok and reactive types. */
function dgReplacer(_key: any, value: any): any {
  if (value instanceof DG.FuncCall) {
    return {
      '#': 'FuncCall',
      id: value.id,
      func: value.func?.nqName,
      inputs: summarizeFuncCallParams([...value.inputParams.values()]),
      outputs: summarizeFuncCallParams([...value.outputParams.values()]),
    };
  }
  if (value instanceof DG.Func)
    return `#Func(${value.nqName})`;
  if (value instanceof DG.DataFrame)
    return summarizeDataFrame(value);
  if (value instanceof Map)
    return Object.fromEntries(value);
  if (value instanceof Set)
    return [...value];
  if (value instanceof BehaviorSubject)
    return value.value;
  if (typeof value === 'function')
    return '#Handler';
  return value;
}

/** Post-process config JSON to reconstruct link spec strings and IO descriptions. */
function humanizeConfig(data: any): any {
  if (data == null || typeof data !== 'object') return data;
  if (Array.isArray(data)) return data.map(humanizeConfig);

  const result: Record<string, any> = {};
  for (const [k, v] of Object.entries(data)) {
    if ((k === 'from' || k === 'to' || k === 'not' || k === 'base' || k === 'actions') && isLinkIOParsedArray(v)) {
      result[k] = linkIOArrayToStrings(v);
    } else if (k === 'io' && Array.isArray(v) && v.length > 0 && v[0]?.id && v[0]?.direction) {
      result[k] = v.map((item: any) => {
        const nullable = item.nullable ? '?' : '';
        return `${item.direction} ${item.id}${nullable}: ${item.type}`;
      });
    } else {
      result[k] = humanizeConfig(v);
    }
  }
  return result;
}

// ---- Tab-specific data builders ----

function toJSON(x: any): any {
  return x ? JSON.parse(JSON.stringify(x, dgReplacer)) : {};
}

/** Build selected step summary from tree state. */
function buildSelectedStepData(
  uuid: string | undefined,
  treeState: PipelineState | undefined,
  stepStates: StepStates | undefined,
): any {
  if (!uuid || !treeState || !stepStates) return undefined;

  const findNode = (state: PipelineState, targetUuid: string): PipelineState | undefined => {
    if (state.uuid === targetUuid) return state;
    if (!isFuncCallState(state)) {
      for (const step of state.steps) {
        const found = findNode(step, targetUuid);
        if (found) return found;
      }
    }
    return undefined;
  };

  const node = findNode(treeState, uuid);
  if (!node) return undefined;

  const data: Record<string, any> = {
    uuid, configId: node.configId, friendlyName: node.friendlyName,
    type: node.type, isReadonly: node.isReadonly,
  };

  if (isFuncCallState(node) && node.funcCall) {
    data.inputs = summarizeFuncCallParams([...node.funcCall.inputParams.values()]);
    data.outputs = summarizeFuncCallParams([...node.funcCall.outputParams.values()]);
  }

  if (stepStates.calls[uuid]) data.callState = stepStates.calls[uuid];
  if (stepStates.validations[uuid]) data.validations = stepStates.validations[uuid];
  if (stepStates.consistency[uuid]) data.consistency = stepStates.consistency[uuid];
  if (stepStates.descriptions[uuid]) data.descriptions = stepStates.descriptions[uuid];
  if (stepStates.meta[uuid]) {
    data.meta = {};
    for (const [k, subj] of Object.entries(stepStates.meta[uuid]!))
      data.meta[k] = subj?.value;
  }

  return data;
}

/** Build human-readable link entries by merging config + runtime data. */
function buildLinksData(links: LinksData[], config: PipelineConfigurationProcessed | undefined) {
  const configLinksMap = new Map<string, any>();
  if (config) {
    const collectLinks = (cfg: any) => {
      for (const link of cfg.links ?? [])
        configLinksMap.set(link.id, link);
      for (const step of cfg.steps ?? cfg.stepTypes ?? []) {
        if (step.links) collectLinks(step);
      }
    };
    collectLinks(config);
  }

  return links.map((linkData) => {
    const spec = configLinksMap.get(linkData.id);
    const matchSpec = linkData.matchInfo?.spec;

    const entry: Record<string, any> = {
      id: linkData.id,
      uuid: linkData.uuid,
      type: spec?.type ?? matchSpec?.type ?? 'data',
      isAction: linkData.isAction,
    };

    const fromParsed = matchSpec?.from ?? spec?.from;
    const toParsed = matchSpec?.to ?? spec?.to;
    if (fromParsed) entry.from = linkIOArrayToStrings(fromParsed);
    if (toParsed) entry.to = linkIOArrayToStrings(toParsed);

    if (spec?.defaultRestrictions) entry.defaultRestrictions = spec.defaultRestrictions;
    if (spec?.handler) entry.handler = '#Handler';
    if (spec?.dataFrameMutations) entry.dataFrameMutations = spec.dataFrameMutations;

    if (linkData.matchInfo?.inputs)
      entry.resolvedInputs = matchedPathsToStrings(linkData.matchInfo.inputs as Record<string, MatchedNodePaths>);
    if (linkData.matchInfo?.outputs)
      entry.resolvedOutputs = matchedPathsToStrings(linkData.matchInfo.outputs as Record<string, MatchedNodePaths>);

    return entry;
  });
}

// ---- Component ----

export const Inspector = Vue.defineComponent({
  name: 'Inspector',
  props: {
    treeState: {
      type: Object as Vue.PropType<PipelineState>,
    },
    links: {
      type: Array as Vue.PropType<LinksData[]>,
    },
    config: {
      type: Object as Vue.PropType<PipelineConfigurationProcessed>,
    },
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
    },
    selectedUuid: {
      type: String,
    },
    stepStates: {
      type: Object as Vue.PropType<StepStates>,
    },
  },
  setup(props) {
    const selectedTab = Vue.ref('Log');
    const linksFilterSelection = Vue.ref<string[]>([]);
    const stepsFilterSelection = Vue.ref<string[]>([]);

    const linksData = Vue.computed(() =>
      buildLinksData(props.links ?? [], props.config));

    // --- Filter options per tab ---

    const linkFilterOptions = Vue.computed<FilterOption[]>(() => {
      const seen = new Set<string>();
      return linksData.value.filter((l) => {
        if (seen.has(l.id)) return false;
        seen.add(l.id);
        return true;
      }).map((l) => ({
        value: l.id,
        label: l.id,
        detail: l.type ?? 'data',
      }));
    });

    const collectSteps = (state: PipelineState): FilterOption[] => {
      const result: FilterOption[] = [];
      const walk = (node: PipelineState) => {
        result.push({
          value: node.uuid,
          label: node.configId,
          detail: node.friendlyName && node.friendlyName !== node.configId ? node.friendlyName : node.type,
        });
        if (!isFuncCallState(node)) {
          for (const step of node.steps)
            walk(step);
        }
      };
      walk(state);
      return result;
    };

    const stepFilterOptions = Vue.computed<FilterOption[]>(() =>
      props.treeState ? collectSteps(props.treeState) : []);

    // --- Filtered data ---

    const handleLinkClicked = (linkId: string) => {
      linksFilterSelection.value = [linkId];
      selectedTab.value = 'Links';
    };

    const filteredLinks = Vue.computed(() => {
      if (!linksFilterSelection.value.length) return linksData.value;
      const sel = new Set(linksFilterSelection.value);
      return linksData.value.filter((l: any) => sel.has(l.id));
    });

    const filterTreeState = (state: PipelineState, uuids: Set<string>): any => {
      if (uuids.has(state.uuid))
        return state;
      if (!isFuncCallState(state)) {
        const filtered = state.steps
          .map((s) => filterTreeState(s, uuids))
          .filter((s) => s != null);
        if (filtered.length > 0)
          return {...state, steps: filtered};
      }
      return undefined;
    };

    const filteredTreeState = Vue.computed(() => {
      if (!props.treeState) return undefined;
      if (!stepsFilterSelection.value.length) return props.treeState;
      return filterTreeState(props.treeState, new Set(stepsFilterSelection.value)) ?? props.treeState;
    });

    const filteredConfig = Vue.computed(() => {
      if (!props.config) return undefined;
      if (!stepsFilterSelection.value.length) return props.config;
      // For config, reuse selected step data when a single step is selected
      if (stepsFilterSelection.value.length === 1)
        return buildSelectedStepData(stepsFilterSelection.value[0], props.treeState, props.stepStates) ?? props.config;
      return props.config;
    });

    const height = 'calc(100% - 20px)';
    const sectionStyle = {height, display: 'flex', flexDirection: 'column' as const};

    const lastVisibleIdx = Vue.ref(0);
    return () => (
      <div style={{userSelect: 'text', overflow: 'hidden'}}>
        <div style={{display: 'flex', gap: '5px', alignItems: 'center', flexWrap: 'wrap'}}>
          <select v-model={selectedTab.value}>
            <option>Log</option>
            <option>Tree State</option>
            <option>Links</option>
            <option>Config</option>
          </select>
          { selectedTab.value === 'Links' &&
            <FilterDropdown
              options={linkFilterOptions.value}
              modelValue={linksFilterSelection.value}
              onUpdate:modelValue={(v: string[]) => linksFilterSelection.value = v}
              placeholder='all links'
            />
          }
          { (selectedTab.value === 'Tree State' || selectedTab.value === 'Config') &&
            <FilterDropdown
              options={stepFilterOptions.value}
              modelValue={stepsFilterSelection.value}
              onUpdate:modelValue={(v: string[]) => stepsFilterSelection.value = v}
              placeholder='all steps'
            />
          }
        </div>
        { selectedTab.value === 'Log' && props.logs &&
          <div style={sectionStyle}>
            <div style={{display: 'flex', flexDirection: 'row'}}>
              <Button onClick={() => lastVisibleIdx.value = 0}>Show All</Button>
              <Button onClick={() => lastVisibleIdx.value = props.logs?.length ?? 0}>Hide Current</Button>
            </div>
            <Logger
              linkFilterOptions={linkFilterOptions.value}
              logs={props.logs.slice(lastVisibleIdx.value)}
              onLinkClicked={handleLinkClicked}
            ></Logger>
          </div>
        }
        { selectedTab.value === 'Tree State' && props.treeState &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={toJSON(filteredTreeState.value)}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Links' && props.links &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={toJSON(filteredLinks.value)}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Config' && props.config &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={humanizeConfig(toJSON(filteredConfig.value))}></VueJsonPretty>
          </div>
        }
      </div>
    );
  },
});
