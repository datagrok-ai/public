import * as Vue from 'vue';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {Tabs, TabItem} from '@datagrok-libraries/webcomponents-vue/src/Tabs/Tabs';

interface Harness {
  root: HTMLElement;
  items: Vue.Ref<TabItem[]>;
  active: Vue.Ref<string>;
  closed: string[];
  added: () => number;
  mounts: Record<string, number>;
  unmounts: Record<string, number>;
  unmount: () => void;
}

const count = (counts: Record<string, number>, name: string) => counts[name] = (counts[name] ?? 0) + 1;

// a pane counts its mounts and keeps local state, to tell a hidden pane from a re-created one
const Pane = Vue.defineComponent({
  props: {
    name: {type: String, required: true},
    mounts: {type: Object as Vue.PropType<Record<string, number>>, required: true},
    unmounts: {type: Object as Vue.PropType<Record<string, number>>, required: true},
  },
  setup(props) {
    const clicks = Vue.ref(0);
    Vue.onMounted(() => count(props.mounts, props.name));
    Vue.onUnmounted(() => count(props.unmounts, props.name));
    return () =>
      Vue.h('button', {class: 'pane-button', onClick: () => clicks.value++}, `${props.name}:${clicks.value}`);
  },
});

function mount(initial: TabItem[], addable = false): Harness {
  const items = Vue.ref(initial);
  const active = Vue.ref(initial[0]?.key ?? '');
  const closed: string[] = [];
  let addCount = 0;
  const mounts: Record<string, number> = {};
  const unmounts: Record<string, number> = {};
  const root = document.createElement('div');
  const app = Vue.createApp({
    render: () => Vue.h(Tabs, {
      'items': items.value,
      'activeKey': active.value,
      'addable': addable,
      'onUpdate:activeKey': (key: string) => active.value = key,
      'onClose': (key: string) => closed.push(key),
      'onAdd': () => addCount++,
    }, {pane: ({item}: {item: TabItem}) => Vue.h(Pane, {name: item.key, mounts, unmounts})}),
  });
  app.config.warnHandler = () => {};
  app.mount(root);
  return {root, items, active, closed, added: () => addCount, mounts, unmounts, unmount: () => app.unmount()};
}

const headers = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLElement>('[role="tab"]'));
const panes = (root: HTMLElement) => Array.from(root.querySelectorAll<HTMLElement>('[role="tabpanel"]'));
const press = (el: HTMLElement, key: string) =>
  el.dispatchEvent(new KeyboardEvent('keydown', {key, bubbles: true, cancelable: true}));
const shownPanes = (root: HTMLElement) => panes(root).filter((pane) => pane.style.display !== 'none')
  .map((pane) => pane.getAttribute('name'));

category('WebComponentsVue: Tabs', () => {
  test('headers carry the platform classes, names and selection', async () => {
    const h = mount([{key: 'a', title: 'First'}, {key: 'b', title: 'Second'}]);
    try {
      expectDeepEqual(headers(h.root).map((el) => el.getAttribute('name')), ['First', 'Second']);
      expect(headers(h.root).every((el) => el.classList.contains('d4-tab-header')), true);
      expectDeepEqual(headers(h.root).map((el) => el.classList.contains('selected')), [true, false]);
      expectDeepEqual(headers(h.root).map((el) => el.getAttribute('aria-selected')), ['true', 'false']);
      expect(h.root.querySelector('.d4-tab-host > .d4-tab-header-stripe > [role="tablist"]') != null, true);
    } finally {
      h.unmount();
    }
  });

  test('a header click shows its pane and keeps the others mounted', async () => {
    const h = mount([{key: 'a', title: 'First'}, {key: 'b', title: 'Second'}]);
    try {
      expectDeepEqual(shownPanes(h.root), ['First']);
      expect(panes(h.root).length, 2);
      const firstButton = panes(h.root)[0].querySelector<HTMLElement>('.pane-button')!;
      firstButton.click();
      await Vue.nextTick();

      headers(h.root)[1].click();
      await Vue.nextTick();
      expect(h.active.value, 'b');
      expectDeepEqual(shownPanes(h.root), ['Second']);

      headers(h.root)[0].click();
      await Vue.nextTick();
      expectDeepEqual(shownPanes(h.root), ['First']);
      expect(panes(h.root)[0].querySelector('.pane-button') === firstButton, true);
      expect(firstButton.textContent, 'a:1');
      expectDeepEqual(h.mounts, {a: 1, b: 1});
    } finally {
      h.unmount();
    }
  });

  test('closing a tab only emits, without switching tabs', async () => {
    const h = mount([{key: 'a', title: 'First'}, {key: 'b', title: 'Second', closable: true}]);
    try {
      expect(headers(h.root)[0].querySelector('[aria-label="Close"]') == null, true);
      headers(h.root)[1].querySelector<HTMLElement>('[aria-label="Close"]')!.click();
      await Vue.nextTick();
      expectDeepEqual(h.closed, ['b']);
      expect(h.active.value, 'a');
    } finally {
      h.unmount();
    }
  });

  test('the + header emits add only when addable', async () => {
    const plain = mount([{key: 'a', title: 'First'}]);
    try {
      expect(plain.root.querySelector('[aria-label="New tab"]') == null, true);
    } finally {
      plain.unmount();
    }
    const h = mount([{key: 'a', title: 'First'}], true);
    try {
      h.root.querySelector<HTMLElement>('[aria-label="New tab"]')!.click();
      await Vue.nextTick();
      expect(h.added(), 1);
      expect(headers(h.root).length, 1);
    } finally {
      h.unmount();
    }
  });

  test('added items mount once and removed ones unmount', async () => {
    const h = mount([{key: 'a', title: 'First'}]);
    try {
      h.items.value = [...h.items.value, {key: 'b', title: 'Second'}];
      await Vue.nextTick();
      expect(panes(h.root).length, 2);
      h.items.value = h.items.value.filter((item) => item.key !== 'a');
      h.active.value = 'b';
      await Vue.nextTick();
      expectDeepEqual(shownPanes(h.root), ['Second']);
      expect(panes(h.root).length, 1);
      expectDeepEqual(h.mounts, {a: 1, b: 1});
      expectDeepEqual(h.unmounts, {a: 1});
    } finally {
      h.unmount();
    }
  });

  test('a pane taller than the host scrolls', async () => {
    const host = document.createElement('div');
    host.style.cssText = 'display: flex; flex-direction: column; height: 200px; position: fixed; top: 0; left: 0;';
    document.body.appendChild(host);
    const app = Vue.createApp({
      render: () => Vue.h(Tabs, {items: [{key: 'a', title: 'First'}], activeKey: 'a', style: {flexShrink: '1'}},
        {pane: () => Vue.h('div', {class: 'tall', style: {height: '1000px', flexShrink: '0'}})}),
    });
    app.config.warnHandler = () => {};
    app.mount(host);
    try {
      const pane = host.querySelector<HTMLElement>('[role="tabpanel"]')!;
      expect(pane.scrollHeight > pane.clientHeight, true);
      pane.scrollTop = pane.scrollHeight;
      const tall = host.querySelector<HTMLElement>('.tall')!.getBoundingClientRect();
      expect(Math.round(tall.bottom) <= Math.round(host.getBoundingClientRect().bottom) + 1, true);
    } finally {
      app.unmount();
      host.remove();
    }
  });

  test('headers work from the keyboard and name their panes', async () => {
    const h = mount([{key: 'a', title: 'First'}, {key: 'b', title: 'Second'}], true);
    const other = mount([{key: 'a', title: 'First'}]);
    try {
      expectDeepEqual(headers(h.root).map((el) => el.getAttribute('tabindex')), ['0', '0']);
      press(headers(h.root)[1], 'Enter');
      await Vue.nextTick();
      expect(h.active.value, 'b');
      press(headers(h.root)[0], ' ');
      await Vue.nextTick();
      expect(h.active.value, 'a');

      headers(h.root).forEach((header, i) => {
        expect(header.getAttribute('aria-controls'), panes(h.root)[i].id);
        expect(panes(h.root)[i].getAttribute('aria-labelledby'), header.id);
      });
      expect(headers(h.root)[0].id !== headers(other.root)[0].id, true);

      const add = h.root.querySelector<HTMLElement>('[aria-label="New tab"]')!;
      expect(add.closest('[role="tablist"]') == null, true);
      expect(add.getAttribute('role'), 'button');
      expect(add.getAttribute('tabindex'), '0');
      press(add, 'Enter');
      expect(h.added(), 1);
    } finally {
      h.unmount();
      other.unmount();
    }
  });
});

