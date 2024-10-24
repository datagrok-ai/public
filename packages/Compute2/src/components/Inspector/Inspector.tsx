import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {PipelineState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';
import {
  BigButton, Button,
} from '@datagrok-libraries/webcomponents-vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {Logger} from '../Logger/Logger';
// it will register a custom element
import 'pretty-json-custom-element';
import {PipelineConfigurationProcessed} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/config-processing-utils';


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
    config: {
      type: Object as Vue.PropType<PipelineConfigurationProcessed>,
    },
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
    },
  },
  setup(props) {
    const selectedTab = Vue.ref('Log');
    const treeStateStr = props.treeState ? JSON.stringify(props.treeState, replacer) : '';
    const configStr = props.config ? JSON.stringify(props.config, replacer) : '';
    // TODO: fix changing between pretty-json elements
    return () => (
      <div>
        <div>
          <select v-model={selectedTab.value}>
            <option>Log</option>
            <option>Tree State</option>
            <option>Config</option>
          </select>
        </div>
        { selectedTab.value === 'Log' && props.logs &&
          <Logger logs={ props.logs }></Logger>
        }
        { selectedTab.value === 'Tree State' && props.treeState &&
          <div>
            <pretty-json expand="5">{treeStateStr}</pretty-json>
          </div>
        }
        { selectedTab.value === 'Config' && props.config &&
          <div>
            <pretty-json expand="5">{configStr}</pretty-json>
          </div>
        }
      </div>
    );
  }
});
