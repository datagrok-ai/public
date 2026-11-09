import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, delay, test} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {BigButton, Button, ComboPopup, IconFA, ValidationIcon} from '@datagrok-libraries/webcomponents-vue';
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

  test('a clickable icon takes focus and presses on Enter and Space', async () => {
    await awaitWebComponents();
    const disabled = Vue.ref(false);
    const clicks: string[] = [];
    let outerKeys = 0;
    const h = mount(() => Vue.h('div', {onKeydown: () => outerKeys++}, [
      Vue.h(IconFA, {name: 'save', isDisabled: disabled.value, onClick: (e: MouseEvent) => clicks.push(e.type)}),
      Vue.h(IconFA, {name: 'lock'}),
    ]));
    const [clickable, plain] = [...h.root.querySelectorAll('dg-icon-fa')] as HTMLElement[];
    // Vue's listeners skip an event from the millisecond they were attached in
    await delay(2);
    try {
      expectDeepEqual(clickable.tabIndex, 0, {prefix: 'clickable'});
      expectDeepEqual(plain.hasAttribute('tabindex'), false, {prefix: 'plain'});
      for (const key of ['Enter', ' ', 'a'])
        clickable.dispatchEvent(new KeyboardEvent('keydown', {key, bubbles: true}));
      expectDeepEqual(clicks, ['click', 'click'], {prefix: 'presses'});
      expectDeepEqual(outerKeys, 1, {prefix: 'only the unhandled key bubbles'});
      disabled.value = true;
      await Vue.nextTick();
      expectDeepEqual(clickable.hasAttribute('tabindex'), false, {prefix: 'disabled'});
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

  test('a disabled big button does not click, says why and stays focusable', async () => {
    await awaitWebComponents();
    const bound = new Map<HTMLElement, unknown>();
    const original = ui.tooltip.bind;
    ui.tooltip.bind = ((el: HTMLElement, value: unknown) => {
      bound.set(el, value);
      return original(el, value as any);
    }) as typeof ui.tooltip.bind;
    const disabled = Vue.ref(true);
    let clicks = 0;
    const h = mount(() => Vue.h(BigButton,
      {isDisabled: disabled.value, disabledReason: 'Fix the inputs', onClick: () => clicks++}, () => 'Run'));
    const button = h.root.querySelector('button')!;
    const tooltipText = () => {
      const value = bound.get(button);
      return typeof value === 'function' ? value() : value;
    };
    const balloons = () => Array.from(document.querySelectorAll('.d4-balloon')).map((el) => el.textContent ?? '');
    try {
      document.querySelectorAll('.d4-balloon').forEach((el) => el.remove());
      expectDeepEqual(button.getAttribute('aria-disabled'), 'true', {prefix: 'aria-disabled'});
      expectDeepEqual(button.tabIndex, 0, {prefix: 'focusable'});
      expectDeepEqual(tooltipText(), 'Fix the inputs', {prefix: 'tooltip'});
      button.click();
      expectDeepEqual(clicks, 0, {prefix: 'no click while disabled'});
      expectDeepEqual(balloons().some((text) => text.includes('Fix the inputs')), true, {prefix: 'warning'});

      disabled.value = false;
      await Vue.nextTick();
      expectDeepEqual(button.hasAttribute('aria-disabled'), false, {prefix: 'enabled aria'});
      expectDeepEqual(tooltipText() ?? null, null, {prefix: 'no tooltip when enabled'});
      button.click();
      expectDeepEqual(clicks, 1, {prefix: 'clicks when enabled'});
    } finally {
      ui.tooltip.bind = original;
      h.dispose();
    }
  });

  test('a button keeps its button classes and passes attributes through', async () => {
    await awaitWebComponents();
    const cls = Vue.ref('a');
    const h = mount(() => Vue.h(Button, {class: cls.value, title: 'Go'}, () => 'Go'));
    const button = () => h.root.querySelector('button')!;
    const classes = () => button().className.split(' ').filter((c) => c).sort();
    try {
      await Vue.nextTick();
      expectDeepEqual(classes(), ['a', 'ui-btn', 'ui-btn-ok'], {prefix: 'mounted'});
      expectDeepEqual(button().title, 'Go', {prefix: 'title'});
      cls.value = 'b';
      await Vue.nextTick();
      expectDeepEqual(classes(), ['b', 'ui-btn', 'ui-btn-ok'], {prefix: 'class changed'});
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
