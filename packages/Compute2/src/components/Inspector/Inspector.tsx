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
import {
  InspectedNode, inspectConfig, LinksInspection, toInspectorJSON,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/inspection';
import {FuncCallStateInfo, ConsistencyInfo} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/StateTreeNodes';
import {ValidationResult} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/common-types';
import {BehaviorSubject} from 'rxjs';
import {FilterDropdown, FilterOption} from './FilterDropdown';

interface StepStates {
  calls: Record<string, FuncCallStateInfo | undefined>;
  validations: Record<string, Record<string, ValidationResult> | undefined>;
  consistency: Record<string, Record<string, ConsistencyInfo> | undefined>;
  meta: Record<string, Record<string, BehaviorSubject<any>> | undefined>;
  descriptions: Record<string, Record<string, string | string[]> | undefined>;
  pipelineValidations: Record<string, ValidationResult | undefined>;
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
    inspectLinks: {
      type: Function as Vue.PropType<() => LinksInspection>,
    },
    inspectNode: {
      type: Function as Vue.PropType<(uuid: string) => InspectedNode | undefined>,
    },
  },
  setup(props) {
    const selectedTab = Vue.ref('Log');
    const linksFilterSelection = Vue.ref<string[]>([]);
    const stepsFilterSelection = Vue.ref<string[]>([]);

    // the driver is read only while the Links tab is shown; links updates after each tree update trigger a re-read
    const linksData = Vue.computed<LinksInspection>(() =>
      selectedTab.value === 'Links' && props.links && props.inspectLinks ?
        props.inspectLinks() :
        {matched: [], notMatched: []});

    // --- Filter options per tab ---

    // the Log tab needs the options too, so matched links come from props.links
    const linkFilterOptions = Vue.computed<FilterOption[]>(() => {
      const seen = new Set<string>();
      const matched = (props.links ?? []).map((l) =>
        ({id: l.id, isAction: l.isAction, type: l.matchInfo.spec.type ?? 'data'}));
      return [...matched, ...linksData.value.notMatched].filter((l) => {
        if (seen.has(l.id)) return false;
        seen.add(l.id);
        return true;
      }).map((l) => ({
        value: l.id,
        label: l.id,
        detail: l.isAction ? `${l.type} (action)` : l.type,
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

    const handleLinkClicked = (linkIds: string[]) => {
      linksFilterSelection.value = [...linkIds];
      selectedTab.value = 'Links';
    };

    const filteredLinks = Vue.computed(() => {
      const {matched, notMatched} = linksData.value;
      if (!linksFilterSelection.value.length) return {matched, notMatched};
      const sel = new Set(linksFilterSelection.value);
      return {matched: matched.filter((l) => sel.has(l.id)), notMatched: notMatched.filter((l) => sel.has(l.id))};
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
      if (stepsFilterSelection.value.length === 1) {
        const uuid = stepsFilterSelection.value[0];
        // the driver read is not reactive, so depend on the step's states here
        const states = props.stepStates;
        void [states?.calls[uuid], states?.validations[uuid], states?.consistency[uuid], states?.meta[uuid],
          states?.descriptions[uuid], states?.pipelineValidations[uuid], props.treeState];
        return props.inspectNode?.(uuid) ?? props.config;
      }
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
            <VueJsonPretty deep={4} showLength={true} data={toInspectorJSON(filteredTreeState.value)}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Links' && props.links &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={filteredLinks.value}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Config' && props.config &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={inspectConfig(filteredConfig.value)}></VueJsonPretty>
          </div>
        }
      </div>
    );
  },
});
