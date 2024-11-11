import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {Button} from '@datagrok-libraries/webcomponents-vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {Logger} from '../Logger/Logger';
import {PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';
import VueJsonPretty from 'vue-json-pretty';
import 'vue-json-pretty/lib/styles.css';
import {LinksData} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/LinksState';

function replacer(key: any, value: any) {
  // Filtering out properties
  if (value instanceof DG.FuncCall) {
    return `#FuncCall`;
  }
  if (value instanceof DG.DataFrame)
    return `#DataFrame`;
  return value;
}

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
  },
  setup(props) {
    const selectedTab = Vue.ref('Log');
    const toJSON = (x: any) => {
      return x ? JSON.parse(JSON.stringify(x, replacer)) : {};
    }

    const lastVisibleIdx = Vue.ref(0);
    return () => (
      <div>
        <div>
          <select v-model={selectedTab.value}>
            <option>Log</option>
            <option>Tree State</option>
            <option>Links</option>
            <option>Config</option>
          </select>
        </div>
        { selectedTab.value === 'Log' && props.logs &&
          <div>
            <div style={{display: 'flex', flexDirection: 'row'}}>
              <Button onClick={() => lastVisibleIdx.value = 0}>Show All</Button>
              <Button onClick={() => lastVisibleIdx.value = props.logs?.length ?? 0}>Hide Current</Button>
            </div>
            <Logger logs={ props.logs.slice(lastVisibleIdx.value) }></Logger>
          </div>
        }
        { selectedTab.value === 'Tree State' && props.treeState &&
          <div>
            <VueJsonPretty deep={3} showLength={true} data={toJSON(props.treeState)}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Links' && props.links &&
          <div>
            <VueJsonPretty deep={3} showLength={true} data={toJSON(props.links)}></VueJsonPretty>
          </div>
        }
        { selectedTab.value === 'Config' && props.config &&
          <div>
            <VueJsonPretty deep={3} showLength={true} data={toJSON(props.config)}></VueJsonPretty>
          </div>
        }
      </div>
    );
  }
});
