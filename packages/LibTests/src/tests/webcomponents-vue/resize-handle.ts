import * as Vue from 'vue';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import {ResizeHandle} from '@datagrok-libraries/webcomponents-vue';

function mount(props: Record<string, any>) {
  const root = document.createElement('div');
  document.body.append(root);
  const size = Vue.ref<number>(props.size ?? 100);
  const updates: number[] = [];
  const ends: number[] = [];
  const app = Vue.createApp({
    render: () => Vue.h(ResizeHandle, {
      'axis': 'x', 'min': 50, 'max': 300, ...props, 'size': size.value,
      'onUpdate:size': (value: number) => {
        updates.push(value);
        size.value = value;
      },
      'onResizeEnd': (value: number) => ends.push(value),
    }),
  });
  app.mount(root);
  const handle = root.firstElementChild as HTMLElement;
  // synthetic pointers have no active pointer to capture
  handle.setPointerCapture = () => {};
  handle.releasePointerCapture = () => {};
  const pointer = (type: string, clientX: number, button = 0) =>
    handle.dispatchEvent(new PointerEvent(type, {clientX, button, pointerId: 1, bubbles: true, cancelable: true}));
  const ghost = () => handle.firstElementChild as HTMLElement | null;
  return {updates, ends, pointer, ghost, dispose: () => {
    app.unmount();
    root.remove();
  }};
}

const nextFrame = () => new Promise((resolve) => requestAnimationFrame(resolve));

category('WebComponentsVue: ResizeHandle', () => {
  test('a drag moves a ghost line and reports the size once on release', async () => {
    const h = mount({});
    try {
      h.pointer('pointerdown', 0);
      h.pointer('pointermove', 40);
      await Vue.nextTick();
      expect(h.ghost() != null, true, 'the ghost line follows the drag');
      expectDeepEqual(h.updates, []);
      h.pointer('pointerup', 40);
      await Vue.nextTick();
      expect(h.ghost(), null);
      expectDeepEqual({updates: h.updates, ends: h.ends}, {updates: [140], ends: [140]});
    } finally {
      h.dispose();
    }
  });

  test('the size is clamped and reverse inverts the drag', async () => {
    for (const [reverse, dx, expected] of [[false, 1000, 300], [false, -1000, 50], [true, 40, 60]] as const) {
      const h = mount({reverse});
      try {
        h.pointer('pointerdown', 0);
        h.pointer('pointermove', dx);
        h.pointer('pointerup', dx);
        expectDeepEqual(h.ends, [expected], {prefix: `reverse ${reverse}, dx ${dx}`});
      } finally {
        h.dispose();
      }
    }
  });

  test('live mode reports once per frame', async () => {
    const h = mount({live: true});
    try {
      h.pointer('pointerdown', 0);
      h.pointer('pointermove', 10);
      h.pointer('pointermove', 20);
      await nextFrame();
      await nextFrame();
      expectDeepEqual(h.updates, [120]);
      expect(h.ghost(), null, 'live mode has no ghost line');
      h.pointer('pointerup', 20);
      expectDeepEqual({updates: h.updates, ends: h.ends}, {updates: [120, 120], ends: [120]});
    } finally {
      h.dispose();
    }
  });

  test('a cancelled pointer ends the drag', async () => {
    const h = mount({});
    try {
      h.pointer('pointerdown', 0);
      h.pointer('pointermove', 40);
      h.pointer('pointercancel', 40);
      await Vue.nextTick();
      expect(h.ghost(), null, 'the ghost line stayed after the pointer was cancelled');
      h.pointer('pointermove', 80);
      await Vue.nextTick();
      expect(h.ghost(), null, 'the handle still follows the pointer');
    } finally {
      h.dispose();
    }
  });

  test('only the primary button starts a drag', async () => {
    const h = mount({});
    try {
      h.pointer('pointerdown', 0, 2);
      h.pointer('pointermove', 40, 2);
      await Vue.nextTick();
      expect(h.ghost(), null, 'a right button press started a drag');
    } finally {
      h.dispose();
    }
  });
});
