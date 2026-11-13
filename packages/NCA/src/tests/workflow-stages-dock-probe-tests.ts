/**
 * nca-workflow-stages S0 — build-entry LIVE probes of the two platform capabilities the stage UI
 * rests on (AC-S0.1, AC-S0.2). Measured on a real `TableView` (rule 16 / evidence class: live
 * measurement), because neither can be settled by reading typings.
 *
 * - AC-S0.1 dock hide/restore: can a docked panel be CLOSED and RE-DOCKED with the SAME element,
 *   its input state and its subscriptions intact? (Stage-scoped widgets = hide/restore of the
 *   existing panels, never an input rebuild — the unit-picker ABANDON lesson.)
 * - AC-S0.2 view-level strip: can a strip sit ABOVE the dock-manager content of a TableView and
 *   survive dock re-layouts and view re-activation?
 *
 * Kept as regression pins: if a platform release changes either behaviour, the stage UI's
 * foundation breaks loudly here rather than silently in the view.
 */

import {category, test, expect, delay} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

function hostView(): DG.TableView {
  const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.INT, 'n', [1, 2, 3])]);
  return grok.shell.addTableView(df);
}

category('Workflow stages: S0 dock probes (live)', () => {
  test('AC-S0.1 a docked panel closes and re-docks with the same element, input state and subscriptions', async () => {
    const view = hostView();
    try {
      const input = ui.input.string('Probe', {value: 'kept'});
      let fired = 0;
      const sub = input.onChanged.subscribe(() => {fired++;});
      const pane = ui.divV([input.root]);
      let node = view.dockManager.dock(pane, 'left', null, 'Probe pane', 0.25);
      await delay(50);
      const w0 = pane.getBoundingClientRect().width;
      view.dockManager.close(node);
      await delay(50);
      expect(pane.isConnected, false, 'closed: detached from the document');
      node = view.dockManager.dock(pane, 'left', null, 'Probe pane', 0.25);
      await delay(50);
      expect(pane.isConnected, true, 're-docked: the SAME element is back in the document');
      expect(input.value, 'kept', 'input value survived the close/re-dock');
      input.value = 'changed';
      await delay(20);
      expect(fired > 0, true, 'the pre-close subscription still fires');
      const w1 = pane.getBoundingClientRect().width;
      expect(Math.abs(w1 - w0) <= Math.max(4, w0 * 0.05), true,
        `dock ratio preserved (width before ${w0.toFixed(0)}px, after ${w1.toFixed(0)}px)`);
      sub.unsubscribe();
      view.dockManager.close(node);
    } finally {
      view.close();
    }
  });

  test('AC-S0.2 a strip above the dock content survives re-layouts and view re-activation', async () => {
    const view = hostView();
    const other = hostView();
    try {
      grok.shell.v = view;
      await delay(50);
      const strip = ui.divText('Setup › Review λz › Report');
      strip.style.flex = '0 0 auto';
      // The view root hosts the dock manager's element; the strip is prepended as a sibling so
      // it takes its own row above it.
      view.root.prepend(strip);
      const panel = ui.divText('docked');
      const node = view.dockManager.dock(panel, 'right', null, 'Probe', 0.3);
      await delay(50);
      const stripRect = strip.getBoundingClientRect();
      const gridRect = view.grid.root.getBoundingClientRect();
      expect(strip.isConnected, true, 'strip in the document');
      expect(stripRect.height > 0, true, `strip has height (${stripRect.height.toFixed(0)}px)`);
      expect(stripRect.bottom <= gridRect.top + 1, true,
        `strip sits above the dock content (strip bottom ${stripRect.bottom.toFixed(0)}, grid top ` +
        `${gridRect.top.toFixed(0)})`);
      view.dockManager.close(node);
      await delay(50);
      expect(strip.isConnected && strip.parentElement === view.root, true, 'survives a dock re-layout');
      grok.shell.v = other;
      await delay(50);
      grok.shell.v = view;
      await delay(50);
      expect(strip.isConnected && strip.parentElement === view.root, true, 'survives view re-activation');
      expect(strip.getBoundingClientRect().bottom <= view.grid.root.getBoundingClientRect().top + 1, true,
        'still above the content after re-activation');
    } finally {
      other.close();
      view.close();
    }
  });
});
