import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {PipelineState, isFuncCallState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Button} from '@datagrok-libraries/webcomponents-vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {Logger} from '../Logger/Logger';
import VueJsonPretty from 'vue-json-pretty';
import 'vue-json-pretty/lib/styles.css';
import {LinksData} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';
import {
  LinksInspection, toInspectorJSON,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/inspection';
import {FilterDropdown, FilterOption} from './FilterDropdown';

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
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
    },
    inspectLinks: {
      type: Function as Vue.PropType<() => LinksInspection>,
    },
    inspectConfig: {
      type: Function as Vue.PropType<() => any>,
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

    // the config only changes when another workflow is loaded, which also replaces the tree state
    const configData = Vue.computed(() =>
      selectedTab.value === 'Config' && props.treeState ? props.inspectConfig?.() : undefined);

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
          { selectedTab.value === 'Tree State' &&
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
        { selectedTab.value === 'Config' && configData.value &&
          <div style={{...sectionStyle, overflow: 'scroll'}}>
            <VueJsonPretty deep={4} showLength={true} data={configData.value}></VueJsonPretty>
          </div>
        }
      </div>
    );
  },
});
