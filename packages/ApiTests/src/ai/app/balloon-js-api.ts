import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow} from '../helpers';

// Each balloon node carries class `d4-balloon` plus its severity class (info/warning/error).
// Auto-dismiss fades after ~5s, so counts are asserted immediately after the synchronous show.
const count = (sel: string = '.d4-balloon'): number => document.querySelectorAll(sel).length;

category('AI: App: Balloon', () => {
  // Clear balloons, run body on a clean slate, always clear again.
  async function withBalloons(body: () => Promise<void> | void): Promise<void> {
    DG.Balloon.closeAll();
    try {
      await body();
    } finally {
      expectNoThrow(() => DG.Balloon.closeAll());
    }
  }

  test('info creates one balloon node with text and info class', async () => {
    await withBalloons(() => {
      expect(count(), 0);
      new DG.Balloon().info('AI balloon info');
      expect(count(), 1);
      const node = document.querySelector('.d4-balloon') as HTMLElement;
      expect(node != null, true);
      expect((node.textContent ?? '').indexOf('AI balloon info') >= 0, true);
      expect(node.classList.contains('info'), true);
    });
  });

  test('warning creates a warning-severity balloon', async () => {
    await withBalloons(() => {
      new DG.Balloon().warning('AI balloon warn');
      expect(count('.d4-balloon.warning'), 1);
      expect(count('.d4-balloon.info'), 0);
    });
  });

  test('error creates an error-severity balloon', async () => {
    await withBalloons(() => {
      new DG.Balloon().error('AI balloon err');
      expect(count('.d4-balloon.error'), 1);
    });
  });

  test('shell.info/warning/error route to Balloon', async () => {
    await withBalloons(() => {
      grok.shell.info('shell-a');
      grok.shell.warning('shell-b');
      grok.shell.error('shell-c');
      expect(count(), 3);
      expect(count('.d4-balloon.info'), 1);
      expect(count('.d4-balloon.warning'), 1);
      expect(count('.d4-balloon.error'), 1);
    });
  });

  test('info accepts an HTMLElement content variant', async () => {
    await withBalloons(() => {
      const el = ui.div('', 'balloon-marker');
      new DG.Balloon().info(el);
      expect(count(), 1);
      const node = document.querySelector('.d4-balloon') as HTMLElement;
      expect(node.contains(el) || node.querySelector('.balloon-marker') != null, true);
    });
  });

  test('multiple infos stack then closeAll clears all', async () => {
    await withBalloons(() => {
      for (let i = 1; i <= 4; i++) {
        new DG.Balloon().info('stack ' + i);
        expect(count(), i);
      }
      DG.Balloon.closeAll();
      expect(count(), 0);
      expect(document.querySelector('.d4-balloon-container') != null, true); // container survives closeAll
    });
  });
}, {owner: 'agolovko@datagrok.ai'});
