import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {KeptAlive} from '../components/KeptAlive/KeptAlive';

type Shown = {uuid: string, cacheKey: string, label: string, visible: boolean};

const providedKey = Symbol('provided');

// the parent is keyed by uuid like a TreeWizard step view, so changing the uuid rebuilds it
function mountHost(first: Shown) {
  const shown = Vue.shallowRef(first);
  const counts = {setups: 0, unmounts: 0};
  const pings: string[] = [];
  let ping: (() => void) | undefined;
  let provided: unknown;

  const Content = Vue.defineComponent({
    props: {label: {type: String, required: true}},
    emits: {ping: () => true},
    setup(props, {emit}) {
      counts.setups++;
      ping = () => emit('ping');
      provided = Vue.inject(providedKey);
      Vue.onUnmounted(() => counts.unmounts++);
      return () => Vue.h('span', props.label);
    },
  });

  const Parent = Vue.defineComponent({
    props: {
      uuid: {type: String, required: true},
      cacheKey: {type: String, required: true},
      label: {type: String, required: true},
      visible: {type: Boolean, required: true},
    },
    setup(props) {
      return () => Vue.h('div', props.visible ? [Vue.h(KeptAlive, null, {
        default: () => [Vue.h(Content, {
          key: props.cacheKey, label: props.label, onPing: () => pings.push(props.uuid),
        })],
      })] : []);
    },
  });

  const app = Vue.createApp({setup: () => () => Vue.h(Parent, {key: shown.value.uuid, ...shown.value})});
  app.provide(providedKey, 'from the app');
  const root = ui.div();
  app.mount(root);
  return {
    root, counts, pings,
    ping: () => ping!(),
    provided: () => provided,
    show: async (next: Partial<Shown>) => {
      shown.value = {...shown.value, ...next};
      await Vue.nextTick();
    },
    unmount: () => app.unmount(),
  };
}

category('KeptAlive', () => {
  test('Keeps the instance when the parent is rebuilt with the same key', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    try {
      await host.show({uuid: '2', label: 'second'});
      expect(host.counts.setups, 1);
      expect(host.counts.unmounts, 0);
      expect(host.root.textContent, 'second');
    } finally {
      host.unmount();
    }
  });

  test('Updates the kept instance in place', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    try {
      await host.show({label: 'second'});
      expect(host.counts.setups, 1);
      expect(host.root.textContent, 'second');
    } finally {
      host.unmount();
    }
  });

  test('Sends events to the current parent', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    try {
      await host.show({uuid: '2'});
      host.ping();
      expectDeepEqual(host.pings, ['2']);
    } finally {
      host.unmount();
    }
  });

  test('Gives the kept instance what the app provides', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    try {
      expect(host.provided(), 'from the app');
    } finally {
      host.unmount();
    }
  });

  test('Replaces the instance for another key', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:First', label: 'first', visible: true});
    try {
      await host.show({uuid: '2', cacheKey: 'Pkg:Second', label: 'second'});
      expect(host.counts.setups, 2);
      expect(host.counts.unmounts, 1);
      expect(host.root.textContent, 'second');

      await host.show({cacheKey: 'Pkg:Third', label: 'third'});
      expect(host.counts.setups, 3);
      expect(host.counts.unmounts, 2);
      expect(host.root.textContent, 'third');
    } finally {
      host.unmount();
    }
  });

  test('Restores the instance when shown again with the same key', async () => {
    const host = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    try {
      await host.show({visible: false});
      expect(host.root.textContent, '');
      expect(host.counts.unmounts, 0);

      await host.show({visible: true, label: 'again'});
      expect(host.counts.setups, 1);
      expect(host.root.textContent, 'again');
    } finally {
      host.unmount();
    }
  });

  test('Destroys the kept instance with the app', async () => {
    const shownHost = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    shownHost.unmount();
    expect(shownHost.counts.unmounts, 1);

    const hiddenHost = mountHost({uuid: '1', cacheKey: 'Pkg:Step', label: 'first', visible: true});
    await hiddenHost.show({visible: false});
    hiddenHost.unmount();
    expect(hiddenHost.counts.unmounts, 1);
  });
});
