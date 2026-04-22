import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {LogItem} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/data/Logger';
import {formatNodePath, formatMutationPath} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/utils';
import {FilterDropdown, FilterOption} from '../Inspector/FilterDropdown';


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
    linkFilterOptions: {
      type: Array as Vue.PropType<FilterOption[]>,
      required: true,
    },
  },
  emits: {
    'linkClicked': (_linkId: string) => true,
  },
  setup(props, {emit}) {
    const eventTypeOptions: FilterOption[] = [
      {value: 'treeUpdateStarted', label: 'treeUpdateStarted', detail: 'tree'},
      {value: 'treeUpdateFinished', label: 'treeUpdateFinished', detail: 'tree'},
      {value: 'treeUpdateMutation', label: 'treeUpdateMutation', detail: 'tree'},
      {value: 'linkAdded', label: 'linkAdded', detail: 'link'},
      {value: 'linkRemoved', label: 'linkRemoved', detail: 'link'},
      {value: 'linkRunStarted', label: 'linkRunStarted', detail: 'link'},
      {value: 'linkRunFinished', label: 'linkRunFinished', detail: 'link'},
      {value: 'actionAdded', label: 'actionAdded', detail: 'action'},
      {value: 'actionRemoved', label: 'actionRemoved', detail: 'action'},
      {value: 'error', label: 'error', detail: 'error'},
    ];

    const eventsFilter = Vue.ref<string[]>([]);
    const linksFilter = Vue.ref<string[]>([]);

    const linkStyle = {
      cursor: 'pointer',
      color: 'var(--blue-1)',
      textDecoration: 'underline',
    };

    const isLinkLogItem = (item: LogItem): item is LogItem & {linkUUID: string; prefix: any; basePath?: any; id: string; isDefaultValidator?: boolean} =>
      item.type === 'linkAdded' || item.type === 'linkRemoved' ||
      item.type === 'linkRunStarted' || item.type === 'linkRunFinished' ||
      item.type === 'actionAdded' || item.type === 'actionRemoved';

    return () => {
      const items = props.logs.filter((item) => {
        if (isLinkLogItem(item) && item.isDefaultValidator)
          return false;
        if (eventsFilter.value.length)
          return eventsFilter.value.includes(item.type);
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
            <div style={{marginRight: '2px'}}>Events:</div>
            <FilterDropdown
              options={eventTypeOptions}
              modelValue={eventsFilter.value}
              onUpdate:modelValue={(v: string[]) => eventsFilter.value = v}
              placeholder='all'
            />
            <div style={{marginRight: '2px', marginLeft: '5px'}}>Links:</div>
            <FilterDropdown
              options={props.linkFilterOptions}
              modelValue={linksFilter.value}
              onUpdate:modelValue={(v: string[]) => linksFilter.value = v}
            />
          </div>
          <div style={{overflow: 'scroll', display: 'grid', gridTemplateColumns: 'auto auto auto', columnGap: '10px', rowGap: '5px'}}>
            {logs}
          </div>
        </>
      );
    };
  },
});
