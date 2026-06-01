import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, withTableView, until, expectNoThrow} from '../helpers';

// DG.JsViewer — public/js-api/src/viewer.ts:379 (members 406-427); Dart host:
// core/client/d4/lib/src/viewers/js_viewer/js_viewer_host_core.dart (scenario: js-viewer)
// JsViewer custom-subclass surface: addRowSourceAndFormula() registering the
// rowSource + filter props in the constructor, and the four lifecycle hooks the
// host fires — onFrameAttached(df) (which chains onTableAttached()),
// onTableAttached(), sourceRowsChanged() (which reads the Dart filter then chains
// onSourceRowsChanged()), and onSourceRowsChanged(). A custom subclass is
// registered via grok.shell.registerViewer(name, desc, factory) and attached with
// tv.addViewer(name) — but addViewer returns the HOST wrapper, so the factory
// closure captures the real subclass instance to assert call counts against.
category('AI: Viewers: JsViewer', () => {
  let counter = 0;
  const uniqueName = () => `Ai-Test-JsViewer-${counter++}`;

  // Inline subclass that records how many times each lifecycle hook fires.
  class TestJsViewer extends DG.JsViewer {
    frameAttachedCount = 0;
    tableAttachedCount = 0;
    sourceRowsChangedCount = 0;
    onSourceRowsChangedCount = 0;

    constructor() {
      super();
      this.addRowSourceAndFormula();
    }

    onFrameAttached(dataFrame: DG.DataFrame): void {
      this.frameAttachedCount++;
      super.onFrameAttached(dataFrame);
    }

    onTableAttached(): void {
      this.tableAttachedCount++;
      super.onTableAttached();
    }

    sourceRowsChanged(): void {
      this.sourceRowsChangedCount++;
      super.sourceRowsChanged();
    }

    onSourceRowsChanged(): void {
      this.onSourceRowsChangedCount++;
      super.onSourceRowsChanged();
    }
  }

  test('addRowSourceAndFormula registers rowSource + filter props on construction', async () => {
    let v: TestJsViewer | undefined;
    expectNoThrow(() => v = new TestJsViewer());
    expect(v!.rowSource, 'Filtered');
    expect(v!.formulaFilter, '');
    const propNames = v!.getProperties().map((p) => p.name);
    expect(propNames.includes('rowSource'), true);
    expect(propNames.includes('filter'), true);
  });

  test('attach fires onFrameAttached + onTableAttached', async () => {
    let captured: TestJsViewer | undefined;
    const name = uniqueName();
    grok.shell.registerViewer(name, 'Ai test JsViewer', () => {
      captured = new TestJsViewer();
      return captured;
    });
    await withTableView(demog(20), async (tv) => {
      tv.addViewer(name);
      await until(() => captured != null && captured!.frameAttachedCount >= 1 && captured!.tableAttachedCount >= 1);
      expect(captured != null, true);
      expect(captured!.frameAttachedCount >= 1, true);
      expect(captured!.tableAttachedCount >= 1, true);
    });
  });

  test('attach fires sourceRowsChanged + onSourceRowsChanged', async () => {
    let captured: TestJsViewer | undefined;
    const name = uniqueName();
    grok.shell.registerViewer(name, 'Ai test JsViewer', () => {
      captured = new TestJsViewer();
      return captured;
    });
    await withTableView(demog(20), async (tv) => {
      tv.addViewer(name);
      await until(() => captured != null && captured!.sourceRowsChangedCount >= 1 && captured!.onSourceRowsChangedCount >= 1);
      expect(captured != null, true);
      expect(captured!.sourceRowsChangedCount >= 1, true);
      expect(captured!.onSourceRowsChangedCount >= 1, true);
    });
  });

  test('post-attach rowSource change re-fires sourceRowsChanged/onSourceRowsChanged', async () => {
    let captured: TestJsViewer | undefined;
    const name = uniqueName();
    grok.shell.registerViewer(name, 'Ai test JsViewer', () => {
      captured = new TestJsViewer();
      return captured;
    });
    await withTableView(demog(20), async (tv) => {
      const host = tv.addViewer(name);
      await until(() => captured != null && captured!.sourceRowsChangedCount >= 1);
      const before = captured!.sourceRowsChangedCount;
      const beforeHook = captured!.onSourceRowsChangedCount;
      // Toggle the row source so the host re-pushes the filtered rows.
      host.setOptions({rowSource: 'Selected'});
      tv.dataFrame.selection.set(0, true);
      tv.dataFrame.selection.set(1, true);
      await until(() => captured!.sourceRowsChangedCount > before && captured!.onSourceRowsChangedCount > beforeHook);
      expect(captured!.sourceRowsChangedCount >= before, true);
      expect(captured!.onSourceRowsChangedCount >= beforeHook, true);
    });
  });

  test('direct-call smoke touches all five members deterministically', async () => {
    // Attach first so the viewer holds a real DataFrame — sourceRowsChanged()
    // reads the Dart combinedFilter, which needs an attached frame.
    let captured: TestJsViewer | undefined;
    const name = uniqueName();
    grok.shell.registerViewer(name, 'Ai test JsViewer', () => {
      captured = new TestJsViewer();
      return captured;
    });
    await withTableView(demog(10), async (tv) => {
      tv.addViewer(name);
      await until(() => captured != null);
      expect(captured != null, true);
      const v = captured!;
      const dfr = tv.dataFrame;
      expectNoThrow(() => {
        v.onFrameAttached(dfr);
        v.onTableAttached();
        v.sourceRowsChanged();
        v.onSourceRowsChanged();
        v.addRowSourceAndFormula();
      });
      // onFrameAttached chains onTableAttached, so both counters advance at least once.
      expect(v.frameAttachedCount >= 1, true);
      expect(v.tableAttachedCount >= 1, true);
      expect(v.sourceRowsChangedCount >= 1, true);
      expect(v.onSourceRowsChangedCount >= 1, true);
    });
  });
});
