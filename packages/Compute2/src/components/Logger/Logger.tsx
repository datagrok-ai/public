import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';


export const Logger = Vue.defineComponent({
  name: 'Logger',
  props: {
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
      required: true,
    },
  },
  setup(props) {
    return () => {
      const logs = props.logs.map((item) => {
        if (item.type === 'treeUpdateStarted' || item.type === 'treeUpdateFinished') {
          const dateString = item.timestamp.toISOString();
          return (
            <div key={item.uuid} style={{display: 'flex', borderBottom: '1px solid Gainsboro', marginTop: '5px'}}>
              <span>
                { dateString } | { item.type }
              </span>
            </div>
          );
        } else if (item.type === 'treeUpdateMutation') {
          const dateString = item.timestamp.toISOString();
          const pathString = item.mutationRootPath?.length ? [...item.mutationRootPath.map(({idx}) => idx)].join('/') : '/';
          return (
            <div key={item.uuid} style={{display: 'flex', borderBottom: '1px solid Gainsboro', marginTop: '5px'}}>
              <span>
                { dateString } | { item.type } | mutationRootPath: { pathString } | addIdx: { String(item.addIdx) } | removeIdx: { String(item.removeIdx) }
              </span>
            </div>
          );
        } else {
          const dateString = item.timestamp.toISOString();
          const pathString = [...item.prefix.map(({idx}) => idx), item.id].join('/');
          return (
            <div key={item.uuid} style={{display: 'flex', borderBottom: '1px solid Gainsboro', marginTop: '5px'}}>
              <span>
                { dateString } | { item.type } | { pathString } | { item.linkUUID }
              </span>
            </div>
          );
        }
      });
      return <div>{ logs }</div>;
    };
  },
});
