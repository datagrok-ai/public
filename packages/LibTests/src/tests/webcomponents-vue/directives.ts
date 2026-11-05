import * as ui from 'datagrok-api/ui';
import * as Vue from 'vue';
import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {ifOverlapping} from '@datagrok-libraries/webcomponents-vue/src/directives/ifOverlapping';
import {wheelGuard} from '@datagrok-libraries/webcomponents-vue/src/directives/wheelGuard';
import {tooltip} from '@datagrok-libraries/webcomponents-vue/src/directives/tooltip';

const HINT = 'Ctrl + scroll to zoom';

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

const wheel = (el: HTMLElement, ctrlKey = false) => {
  const ev = new WheelEvent('wheel', {deltaY: 100, ctrlKey, bubbles: true, cancelable: true});
  el.dispatchEvent(ev);
  return ev;
};

category('WebComponentsVue: directives', () => {
  test('ifOverlapping shows the loader after the debounce and hides it', async () => {
    const shown = Vue.ref(true);
    const h = mount(() => Vue.withDirectives(Vue.h('div', {class: 'target'}), [[ifOverlapping, shown.value, 'Busy']]));
    try {
      const target = h.root.querySelector('.target')!;
      expect(target.querySelector('.d4-update-shadow'), null, 'the loader is debounced');
      await waitUntil(() => target.querySelector('.d4-update-shadow') != null, 'loader not shown', 1000);
      expect(target.textContent?.includes('Busy'), true);
      shown.value = false;
      await waitUntil(() => target.querySelector('.d4-update-shadow') == null, 'loader not hidden', 1000);
    } finally {
      h.dispose();
    }
  });

  test('ifOverlapping never shows the loader after unmount', async () => {
    const h = mount(() => Vue.withDirectives(Vue.h('div', {class: 'target'}), [[ifOverlapping, true]]));
    const target = h.root.querySelector('.target')!;
    h.dispose();
    await delay(300);
    expect(target.querySelector('.d4-update-shadow'), null);
  });

  test('wheelGuard stops plain wheel and lets ctrl wheel through', async () => {
    const h = mount(() => Vue.withDirectives(Vue.h('div', {class: 'target'}, [Vue.h('div', {class: 'inner'})]),
      [[wheelGuard]]));
    try {
      const inner = h.root.querySelector<HTMLElement>('.inner')!;
      const reached: boolean[] = [];
      inner.addEventListener('wheel', (ev) => reached.push(ev.ctrlKey));
      const plain = wheel(inner);
      const zoom = wheel(inner, true);
      expectDeepEqual(reached, [true], {prefix: 'only ctrl wheel reaches the viewer'});
      expect(plain.defaultPrevented, false, 'plain wheel still scrolls the page');
      expect(zoom.defaultPrevented, true, 'ctrl wheel does not zoom the page');
      const hint = Array.from(h.root.querySelectorAll<HTMLElement>('div'))
        .find((el) => !el.childElementCount && el.textContent === HINT)!;
      expect(hint.style.opacity, '1');
    } finally {
      h.dispose();
    }
  });

  test('the tooltip follows its binding', async () => {
    const bound = new Map<HTMLElement, unknown>();
    const original = ui.tooltip.bind;
    ui.tooltip.bind = ((el: HTMLElement, value: unknown) => {
      bound.set(el, value);
      return original(el, value as any);
    }) as typeof ui.tooltip.bind;
    const text = Vue.ref('first');
    const h = mount(() => Vue.withDirectives(Vue.h('div', {class: 'target'}), [[tooltip, text.value]]));
    try {
      const target = h.root.querySelector<HTMLElement>('.target')!;
      const shown = () => {
        const value = bound.get(target);
        return typeof value === 'function' ? value() : value;
      };
      expect(shown(), 'first');
      text.value = 'second';
      await Vue.nextTick();
      expect(shown(), 'second', 'the tooltip kept its first text');
    } finally {
      ui.tooltip.bind = original;
      h.dispose();
    }
  });
});
