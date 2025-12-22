import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {NodePath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/BaseTree';


export const Logger = Vue.defineComponent({
  name: 'Logger',
  props: {
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
      required: true,
    },
    linkIds: {
      type: Array as  Vue.PropType<string[]>,
      required: true,
    },
  },
  setup(props) {
    const typeFilterOpts = ['all', 'tree', 'links'] as const;
    const treeEvents = new Set(['treeUpdateStarted', 'treeUpdateFinished', 'treeUpdateMutation']);

    const typeFilter = Vue.ref<typeof typeFilterOpts[number]>('all');
    const linksFilter = Vue.ref<string[]>([]);
    const linksFilterOpts = Vue.computed(() => props.linkIds);

    const formatPath = (path?: NodePath, addIdx?: number, removeIdx?: number, id?: string) => {
      let str = path?.length ? [...path.map(({idx, id}) => `[${idx}]${id}`)].join('/') : '';
      const append = (suffix : string) => str += str ? `/${suffix}` : `${suffix}`;
      if (addIdx != null && removeIdx != null)
        append(`[${removeIdx}→${addIdx}]${id}`);
      else if (removeIdx != null)
        append(`-[${removeIdx}]${id}`);
      else if (addIdx != null) {
        append(`+[${addIdx}]${id}`);
      } else if (id) {
        append(id);
      }
      return str;
    }

    return () => {
      const items = props.logs.filter((item) => {
        if (typeFilter.value === 'all')
          return true;
        if (typeFilter.value === 'tree')
          return treeEvents.has(item.type);
        if (typeFilter.value === 'links')
          return !treeEvents.has(item.type);
      }).filter((item) => {
        if (linksFilter.value.length) {
          if (item.type === 'linkAdded' || item.type === 'linkRemoved' || item.type === 'linkRunStarted' || item.type === 'linkRunFinished' || item.type === 'actionAdded' || item.type === 'actionRemoved')
            return linksFilter.value.includes(item.id);
        }
        return true;
      });
      const logs = items.flatMap((item) => {
        if (item.type === 'treeUpdateStarted' || item.type === 'treeUpdateFinished') {
          const dateString = item.timestamp.toISOString();
          return ([
            <div key={item.uuid + '1'}>
              {dateString}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
            </div>,
            <div key={item.uuid + '3'}>
            </div>,
            <div key={item.uuid + '4'}>
            </div>
          ]);
        } else if (item.type === 'treeUpdateMutation') {
          const dateString = item.timestamp.toISOString();
          const pathString = formatPath(item.mutationRootPath, item.addIdx, item.removeIdx, item.id);
          return ([
            <div key={item.uuid + '1'}>
              {dateString}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
              </div>,
            <div key={item.uuid + '3'}>
              {pathString}
            </div>,
            <div key={item.uuid + '4'}>
            </div>
          ]);
        } else {
          const dateString = item.timestamp.toISOString();
          const pathString = formatPath([...item.prefix, ...(item.basePath ?? [])], undefined, undefined, item.id);
          return ([
            <div key={item.uuid + '1'}>
              {dateString}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
              </div>,
            <div key={item.uuid + '3'}>
              {pathString}
            </div>,
            <div key={item.uuid + '4'}>
              {item.linkUUID}
            </div>
          ]);
        }
      });
      return (
        <>
          <div style={{display: 'flex', flexDirection: 'row', marginBottom: '10px'}}>
            <div style={{marginRight: '5px'}}>Events:</div>
            <div style={{marginRight: '20px'}}>
              <select onChange={(ev) => typeFilter.value = (ev.target as any)?.value ?? 'all'} value={typeFilter.value}>
                {typeFilterOpts.map(option => (
                  <option value={option}>{option}</option>
                ))}
              </select>
            </div>
            <div style={{marginRight: '5px'}}>Links:</div>
            <div style={{marginRight: '20px'}}>
              <select multiple onChange={(ev) => linksFilter.value = Array.from((ev.target as any)?.selectedOptions ?? []).map((x: any) => x.value)} value={linksFilter.value}>
                {linksFilterOpts.value.map(option => (
                  <option key={option} value={option}>{option}</option>
                ))}
              </select>
            </div>
            <div> Selected: {linksFilter.value.join(', ')} </div>
          </div>
          <div style={{overflow: 'scroll', display: 'grid', gridTemplateColumns: 'auto auto auto auto', columnGap: '10px', rowGap: '5px'}}>
            {logs}
          </div>
        </>
      );
    };
  },
});
