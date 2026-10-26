import * as Vue from 'vue';
import './Tabs.css';

export interface TabItem {
  key: string;
  title: string;
  closable?: boolean;
}

// Vue.useId counts per app, and each view is its own app: ids must not repeat across views
let lastId = 0;

// Enter and Space press a focused header, like a button
const onPress = (action: () => void) => (e: KeyboardEvent) => {
  if (e.key !== 'Enter' && e.key !== ' ')
    return;
  e.preventDefault();
  action();
};

/**
 * Tab strip styled like ui.tabControl. Controlled: bind the shown tab via `activeKey`/`update:activeKey`
 * (`v-model:activeKey`); `close` and `add` (the + header, when `addable`) leave the items to the parent.
 * Every pane stays mounted and inactive ones are hidden, so their state and viewers survive switching.
 * Headers and panes carry the title as `name`, like the platform's tabs.
 */
export const Tabs = Vue.defineComponent({
  name: 'Tabs',
  props: {
    items: {
      type: Array as Vue.PropType<TabItem[]>,
      required: true,
    },
    activeKey: {
      type: String,
      required: true,
    },
    addable: {
      type: Boolean,
      default: false,
    },
  },
  emits: {
    'update:activeKey': (_key: string) => true,
    'close': (_key: string) => true,
    'add': () => true,
  },
  slots: Object as Vue.SlotsType<{pane: {item: TabItem, active: boolean}}>,
  setup(props, {emit, slots}) {
    const id = `webcomponents-vue-tabs-${++lastId}`;
    const tabId = (index: number) => `${id}-tab-${index}`;
    const paneId = (index: number) => `${id}-pane-${index}`;

    return () => Vue.h('div', {class: 'd4-tab-host webcomponents-vue-tabs'}, [
      Vue.h('div', {class: 'd4-tab-header-stripe'}, [
        Vue.h('div', {class: 'webcomponents-vue-tabs-list', role: 'tablist'}, props.items.map((item, index) => {
          const active = item.key === props.activeKey;
          const select = () => emit('update:activeKey', item.key);
          return Vue.h('div', {
            'key': item.key,
            'id': tabId(index),
            'class': active ? 'd4-tab-header selected' : 'd4-tab-header',
            'name': item.title,
            'role': 'tab',
            'tabindex': 0,
            'aria-selected': active,
            'aria-controls': paneId(index),
            'onClick': select,
            'onKeydown': onPress(select),
          }, [
            Vue.h('span', item.title),
            item.closable &&
              Vue.h('span', {class: 'webcomponents-vue-tabs-close', onClick: (e: MouseEvent) => e.stopPropagation()}, [
                Vue.h('dg-icon-fa', {
                  'name': 'times', 'tooltip': 'Close', 'role': 'button', 'aria-label': 'Close',
                  'onClick': () => emit('close', item.key),
                }),
              ]),
          ]);
        })),
        props.addable &&
          Vue.h('div', {
            'class': 'd4-tab-header',
            'role': 'button',
            'tabindex': 0,
            'aria-label': 'New tab',
            'onClick': () => emit('add'),
            'onKeydown': onPress(() => emit('add')),
          }, [Vue.h('dg-icon-fa', {name: 'plus', tooltip: 'New tab'})]),
      ]),
      Vue.h('div', {class: 'd4-tab-content'}, props.items.map((item, index) => {
        const active = item.key === props.activeKey;
        return Vue.h('div', {
          'key': item.key,
          'id': paneId(index),
          'class': 'webcomponents-vue-tabs-pane',
          'name': item.title,
          'role': 'tabpanel',
          'aria-labelledby': tabId(index),
          'style': active ? undefined : {display: 'none'},
        }, slots.pane?.({item, active}));
      })),
    ]);
  },
});
