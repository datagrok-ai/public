import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {BigButton, ComboPopup, IconFA, ValidationIcon} from '@datagrok-libraries/webcomponents-vue';
import {awaitWebComponents} from '../../test-utils';
import {consistency, findLink, openPopover} from '../webcomponents/fixtures';

function mount(render: () => Vue.VNode) {
  const root = ui.div();
  document.body.append(root);
  const app = Vue.createApp({render});
  app.mount(root);
  return {root, dispose: () => {
    app.unmount();
    root.remove();
  }};
}

const naming = (el: Element | null) =>
  ({role: el?.getAttribute('role') ?? null, label: el?.getAttribute('aria-label') ?? null});

category('WebComponentsVue: element wrappers', () => {
  test('an icon is named by its tooltip, as a button only when clickable', async () => {
    await awaitWebComponents();
    const cases: [Record<string, any>, ReturnType<typeof naming>][] = [
      [{tooltip: 'Save', onClick: () => {}}, {role: 'button', label: 'Save'}],
      [{tooltip: 'Locked'}, {role: 'img', label: 'Locked'}],
      [{onClick: () => {}}, {role: null, label: null}],
    ];
    for (const [props, expected] of cases) {
      const h = mount(() => Vue.h(IconFA, {name: 'save', ...props}));
      try {
        expectDeepEqual(naming(h.root.querySelector('dg-icon-fa')), expected,
          {prefix: JSON.stringify(Object.keys(props))});
      } finally {
        h.dispose();
      }
    }
  });

  test('a disabled icon stays disabled when its class changes', async () => {
    await awaitWebComponents();
    const cls = Vue.ref('a');
    const disabled = Vue.ref(true);
    const h = mount(() => Vue.h(IconFA, {name: 'save', class: cls.value, isDisabled: disabled.value}));
    const classes = () => h.root.querySelector('dg-icon-fa')!.className;
    try {
      await Vue.nextTick();
      expectDeepEqual(classes().split(' ').sort(), ['a', 'd4-disabled'], {prefix: 'mounted'});
      cls.value = 'b';
      await Vue.nextTick();
      expectDeepEqual(classes().split(' ').sort(), ['b', 'd4-disabled'], {prefix: 'class changed'});
      disabled.value = false;
      await Vue.nextTick();
      expectDeepEqual(classes(), 'b', {prefix: 'enabled'});
    } finally {
      h.dispose();
    }
  });

  test('a disabled big button keeps its button classes when its class changes', async () => {
    await awaitWebComponents();
    const cls = Vue.ref('a');
    const disabled = Vue.ref(true);
    const h = mount(() => Vue.h(BigButton, {class: cls.value, isDisabled: disabled.value}, () => 'Run'));
    const classes = () => h.root.querySelector('button')!.className.split(' ').filter((c) => c).sort();
    const base = ['ui-btn', 'ui-btn-ok', 'ui-btn-raised'];
    try {
      await Vue.nextTick();
      expectDeepEqual(classes(), ['a', 'd4-disabled', ...base], {prefix: 'mounted'});
      cls.value = 'b';
      await Vue.nextTick();
      expectDeepEqual(classes(), ['b', 'd4-disabled', ...base], {prefix: 'class changed'});
      disabled.value = false;
      await Vue.nextTick();
      expectDeepEqual(classes(), ['b', ...base], {prefix: 'enabled'});
    } finally {
      h.dispose();
    }
  });

  test('a combo popup is named by its tooltip', async () => {
    await awaitWebComponents();
    const h = mount(() => Vue.h(ComboPopup, {caption: 'Pick', items: ['a'], tooltip: 'Choose one'}));
    try {
      expectDeepEqual(naming(h.root.querySelector('dg-combo-popup')), {role: 'button', label: 'Choose one'});
    } finally {
      h.dispose();
    }
  });

  test('the validation icon forwards the dataframe diff request', async () => {
    await awaitWebComponents();
    const requests: unknown[] = [];
    const h = mount(() => Vue.h(ValidationIcon, {
      validationStatus: {consistency: consistency('restricted', true, 'computed')},
      isScalar: false, isDataFrame: true,
      onShowDataFrameDiff: (value: unknown) => requests.push(value),
    }));
    try {
      findLink(await openPopover(h.root), 'Show dataframes diff')!.click();
      expectDeepEqual(requests, ['computed'], {prefix: 'the diff link does nothing through the Vue wrapper'});
    } finally {
      h.dispose();
    }
  });
});
