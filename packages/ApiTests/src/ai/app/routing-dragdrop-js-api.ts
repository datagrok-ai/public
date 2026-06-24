import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {expectNoThrow} from '../helpers';

// Routing + drag-drop public surface: grok.shell.route (shell.ts:192), ui.makeDraggable /
// ui.makeDroppable (ui.ts:723,776), and DG.DragDropArgs.handled (events.ts:303-336).
category('AI: App: Routing and DragDrop JS API', () => {
  test('shell.route returns a View', async () => {
    const prior = grok.shell.v;
    let routed: any = null;
    try {
      routed = grok.shell.route('/');
      expect(routed instanceof DG.View, true);
    } finally {
      // route() may navigate; restore the previously-current view if it still exists.
      if (prior != null && grok.shell.v !== prior)
        expectNoThrow(() => grok.shell.v = prior);
      if (routed != null && routed !== prior && routed.close != null)
        expectNoThrow(() => routed.close());
    }
  });

  test('makeDraggable returns the same element', async () => {
    const el = ui.div('drag-me');
    try {
      const ret = ui.makeDraggable(el);
      expect(ret === el, true);
    } finally {
      el.remove();
    }
  });

  test('makeDraggable with getDragObject wiring does not throw', async () => {
    const el = ui.div('drag-payload');
    try {
      const payload = {kind: 'demo'};
      expectNoThrow(() => ui.makeDraggable(el, {
        getDragObject: () => payload,
        getDragCaption: () => 'Dragging demo',
        dragObjectType: 'demo',
        check: () => true,
        allowCopy: () => true,
      }));
    } finally {
      el.remove();
    }
  });

  test('makeDroppable with acceptDrag + acceptDrop + doDrop wiring does not throw', async () => {
    const el = ui.div('drop-zone');
    try {
      let dropped = false;
      // makeDroppable returns void; we assert the wiring registers without throwing.
      expectNoThrow(() => ui.makeDroppable(el, {
        acceptDrag: (args) => {
          // args is a DG.DragDropArgs; veto by setting handled.
          args.handled = true;
          return false;
        },
        acceptDrop: (o) => o instanceof DG.Column,
        doDrop: (_args) => {dropped = true;},
        dropSuggestion: 'Drop a numeric column',
      }));
      // No synthetic drag dispatched, so the handler must not have fired yet.
      expect(dropped, false);
    } finally {
      el.remove();
    }
  });

  test('DragDropArgs exposes a read/write handled accessor', async () => {
    // A live DragDropArgs only exists during a real drag session; verify the public
    // accessor surface on the prototype instead of fabricating an invalid handle.
    expect(typeof DG.DragDropArgs, 'function');
    const desc = Object.getOwnPropertyDescriptor(DG.DragDropArgs.prototype, 'handled');
    expect(desc != null, true);
    expect(typeof desc!.get, 'function');
    expect(typeof desc!.set, 'function');
  });
}, {owner: 'agolovko@datagrok.ai'});
