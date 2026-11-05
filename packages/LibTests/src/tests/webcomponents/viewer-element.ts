import * as DG from 'datagrok-api/dg';
import {category, test, expect, waitUntil, delay} from '@datagrok-libraries/test/src/test';
import {expectDeepEqual} from '@datagrok-libraries/utils/src/expect';
import type {Viewer} from '@datagrok-libraries/webcomponents/src/Viewer/Viewer';
import {awaitWebComponents} from '../../test-utils';

const makeDf = () => DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'x', [1, 2, 3])]);

const EVENTS = ['viewer-changed', 'viewer-data-frame-changed'];

async function mountViewer() {
  await awaitWebComponents();
  const el = document.createElement('dg-viewer') as Viewer;
  el.style.cssText = 'display:block;width:400px;height:300px';
  const events: string[] = [];
  for (const name of EVENTS)
    el.addEventListener(name, () => events.push(name));
  document.body.append(el);
  const count = (name: string) => events.filter((e) => e === name).length;
  const dispose = () => {
    el.destroy();
    el.remove();
  };
  return {el, events, count, dispose};
}

category('WebComponents: Viewer element', () => {
  test('creation is not announced as a frame swap, an in-place swap once', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      const df1 = makeDf();
      el.type = DG.VIEWER.GRID;
      el.dataFrame = df1;
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      await delay(100);
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 1, swapped: 0});
      const viewer = el.viewer;
      const df2 = makeDf();
      el.dataFrame = df2;
      await waitUntil(() => el.viewer?.dataFrame?.dart === df2.dart, 'frame not swapped', 5000);
      await delay(100);
      expect(el.viewer === viewer, true, 'a frame swap keeps the viewer');
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 1, swapped: 1});
    } finally {
      dispose();
    }
  });

  test('an options change is not announced as a frame swap', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      el.options = {showRowHeader: false};
      await delay(300);
      expect(count('viewer-data-frame-changed'), 0);
    } finally {
      dispose();
    }
  });

  test('a frame set on the viewer is announced once, outside the Dart event, and kept', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      const df1 = makeDf();
      el.type = DG.VIEWER.GRID;
      el.dataFrame = df1;
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      await delay(100);
      let settingFrame = false;
      let insideSetter = false;
      el.addEventListener('viewer-data-frame-changed', () => insideSetter ||= settingFrame);
      const df2 = makeDf();
      settingFrame = true;
      el.viewer!.dataFrame = df2;
      settingFrame = false;
      await delay(300);
      expect(count('viewer-data-frame-changed'), 1);
      expect(insideSetter, false, 'announced inside the Dart event');
      el.options = {showRowHeader: false};
      await delay(300);
      expect(el.viewer!.dataFrame.dart === df2.dart, true, 'an options change restored the old frame');
    } finally {
      dispose();
    }
  });

  test('a type change keeping the frame announces only the viewer', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      el.type = DG.VIEWER.SCATTER_PLOT;
      await waitUntil(() => el.viewer?.type === DG.VIEWER.SCATTER_PLOT, 'viewer not replaced', 5000);
      await delay(300);
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 2, swapped: 0});
    } finally {
      dispose();
    }
  });

  test('a frame without a type announces nothing until the viewer is created', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      el.dataFrame = makeDf();
      await delay(300);
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 0, swapped: 0}, {prefix: 'without a type'});
      el.type = DG.VIEWER.GRID;
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      await delay(100);
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 1, swapped: 0});
    } finally {
      dispose();
    }
  });

  test('clearing the frame announces only the viewer', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      el.dataFrame = undefined;
      await waitUntil(() => !el.viewer, 'viewer not detached', 5000);
      await delay(100);
      expectDeepEqual({changed: count('viewer-changed'), swapped: count('viewer-data-frame-changed')},
        {changed: 2, swapped: 0});
    } finally {
      dispose();
    }
  });

  test('destroy empties the element and ignores later properties', async () => {
    const {el, count, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      el.destroy();
      expect(el.childElementCount, 0);
      el.type = DG.VIEWER.SCATTER_PLOT;
      el.dataFrame = makeDf();
      await delay(300);
      expect(el.childElementCount, 0);
      expect(count('viewer-changed'), 1);
    } finally {
      dispose();
    }
  });

  test('an unknown viewer type does not break the element', async () => {
    const {el, dispose} = await mountViewer();
    try {
      el.type = 'No such viewer';
      el.dataFrame = makeDf();
      await delay(300);
      el.type = DG.VIEWER.GRID;
      await waitUntil(() => el.viewer?.type === DG.VIEWER.GRID, 'the element stopped after the failed creation', 3000);
    } finally {
      dispose();
    }
  });

  test('a replaced viewer is detached', async () => {
    const {el, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      let detached = false;
      el.viewer!.onDetached.subscribe(() => detached = true);
      el.type = DG.VIEWER.SCATTER_PLOT;
      await waitUntil(() => el.viewer?.type === DG.VIEWER.SCATTER_PLOT, 'viewer not replaced', 5000);
      await delay(200);
      expect(detached, true, 'the grid replaced by the scatter plot was never detached');
    } finally {
      dispose();
    }
  });

  test('destroy detaches the viewer', async () => {
    const {el, dispose} = await mountViewer();
    try {
      el.type = DG.VIEWER.GRID;
      el.dataFrame = makeDf();
      await waitUntil(() => !!el.viewer, 'viewer not created', 5000);
      let detached = false;
      el.viewer!.onDetached.subscribe(() => detached = true);
      el.destroy();
      await delay(200);
      expect(detached, true, 'the destroyed element left its viewer attached');
    } finally {
      dispose();
    }
  });
});
