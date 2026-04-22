import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {formatNodePath, formatMutationPath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';


function formatTime(d: Date): string {
  const h = String(d.getHours()).padStart(2, '0');
  const m = String(d.getMinutes()).padStart(2, '0');
  const s = String(d.getSeconds()).padStart(2, '0');
  const ms = String(d.getMilliseconds()).padStart(3, '0');
  return `${h}:${m}:${s}.${ms}`;
}

export const Logger = Vue.defineComponent({
  name: 'Logger',
  props: {
    logs: {
      type: Array as Vue.PropType<LogItem[]>,
      required: true,
    },
    linkIds: {
      type: Array as Vue.PropType<string[]>,
      required: true,
    },
  },
  emits: {
    'linkClicked': (_linkId: string) => true,
  },
  setup(props, {emit}) {
    const typeFilterOpts = ['all', 'tree', 'links'] as const;
    const treeEvents = new Set(['treeUpdateStarted', 'treeUpdateFinished', 'treeUpdateMutation']);

    const typeFilter = Vue.ref<typeof typeFilterOpts[number]>('all');
    const linksFilter = Vue.ref<string[]>([]);
    const linksFilterOpts = Vue.computed(() => props.linkIds);
    const showDefaultValidators = Vue.ref(false);

    const linkStyle = {
      cursor: 'pointer',
      color: 'var(--blue-1)',
      textDecoration: 'underline',
    };


    const isLinkLogItem = (item: LogItem): item is LogItem & {linkUUID: string; prefix: any; basePath?: any; id: string; isDefaultValidator?: boolean} =>
      item.type === 'linkAdded' || item.type === 'linkRemoved' ||
      item.type === 'linkRunStarted' || item.type === 'linkRunFinished' ||
      item.type === 'actionAdded' || item.type === 'actionRemoved';

    const checkboxStyle = {display: 'flex', alignItems: 'center', gap: '3px', fontSize: '12px', whiteSpace: 'nowrap' as const};

    return () => {
      const items = props.logs.filter((item) => {
        if (!showDefaultValidators.value && isLinkLogItem(item) && item.isDefaultValidator)
          return false;
        if (typeFilter.value === 'tree')
          return treeEvents.has(item.type);
        if (typeFilter.value === 'links')
          return !treeEvents.has(item.type);
        return true;
      }).filter((item) => {
        if (linksFilter.value.length) {
          if (isLinkLogItem(item))
            return linksFilter.value.includes(item.id);
        }
        return true;
      });
      const logs = items.flatMap((item) => {
        if (item.type === 'treeUpdateStarted' || item.type === 'treeUpdateFinished') {
          return ([
            <div key={item.uuid + '1'}>
              {formatTime(item.timestamp)}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
            </div>,
            <div key={item.uuid + '3'}>
            </div>,
          ]);
        } else if (item.type === 'treeUpdateMutation') {
          const pathString = formatMutationPath(item.mutationRootPath, item.addIdx, item.removeIdx, item.id);
          return ([
            <div key={item.uuid + '1'}>
              {formatTime(item.timestamp)}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
              </div>,
            <div key={item.uuid + '3'}>
              {pathString}
            </div>,
          ]);
        } else if (item.type === 'error') {
          return ([
            <div key={item.uuid + '1'}>
              {formatTime(item.timestamp)}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type} ({item.severity})
              </div>,
            <div key={item.uuid + '3'}>
              {item.context}: {item.message}
            </div>,
          ]);
        } else {
          const segments = [...item.prefix, ...(item.basePath ?? [])].map((s) => ({idx: s.idx, id: s.id}));
          const pathStr = formatNodePath(segments);
          const pathString = pathStr ? `${pathStr}/${item.id}` : item.id;
          const isDefault = item.isDefaultValidator;
          return ([
            <div key={item.uuid + '1'}>
              {formatTime(item.timestamp)}
            </div>,
            <div key={item.uuid + '2'}>
              {item.type}
              </div>,
            <div key={item.uuid + '3'}>
              { isDefault
                ? <span>{pathString}</span>
                : <span style={linkStyle} onClick={() => emit('linkClicked', item.id)}>{pathString}</span>
              }
            </div>,
          ]);
        }
      });
      return (
        <>
          <div style={{display: 'flex', flexDirection: 'row', marginBottom: '10px', alignItems: 'center', flexWrap: 'wrap', gap: '5px'}}>
            <div style={{marginRight: '5px'}}>Events:</div>
            <div style={{marginRight: '10px'}}>
              <select onChange={(ev) => typeFilter.value = (ev.target as any)?.value ?? 'all'} value={typeFilter.value}>
                {typeFilterOpts.map(option => (
                  <option value={option}>{option}</option>
                ))}
              </select>
            </div>
            <div style={{marginRight: '5px'}}>Links:</div>
            <div style={{marginRight: '10px'}}>
              <select multiple onChange={(ev) => linksFilter.value = Array.from((ev.target as any)?.selectedOptions ?? []).map((x: any) => x.value)} value={linksFilter.value}>
                {linksFilterOpts.value.map(option => (
                  <option key={option} value={option}>{option}</option>
                ))}
              </select>
            </div>
            <div style={{marginRight: '10px'}}> Selected: {linksFilter.value.join(', ')} </div>
            <label style={checkboxStyle}>
              <input type='checkbox' v-model={showDefaultValidators.value} />
              Default validators
            </label>
          </div>
          <div style={{overflow: 'scroll', display: 'grid', gridTemplateColumns: 'auto auto auto', columnGap: '10px', rowGap: '5px'}}>
            {logs}
          </div>
        </>
      );
    };
  },
});
