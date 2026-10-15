import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';

category('U2', () => {
  test('Control.forElement', async () => {
    const control = new DG.U2.Control();
    const child = document.createElement('span');
    control.root.append(child);
    expect(DG.U2.Control.forElement(child), control);
    expect(DG.U2.Control.forElement(control.root), control);
    expect(DG.U2.Control.forElement(document.createElement('div')) === undefined, true);

    const inner = new DG.U2.Control();
    control.root.append(inner.root);
    expect(DG.U2.Control.forElement(inner.root), inner);

    const widget = DG.Widget.fromRoot(ui.divText('w'));
    expect(DG.U2.Control.forElement(widget.root), widget);
    widget.detach();
    inner.dispose();
    control.dispose();
  }, {owner: 'askalkin@datagrok.ai'});

  test('getProperties accessor tier', async () => {
    class Sample extends DG.U2.Control {
      readonly value = DG.U2.signal('a');
      private _title = 'untitled';
      readonly fixed = 'construction-only';

      get title(): string {
        return this._title;
      }

      set title(x: string) {
        this._title = x;
        this.root.dataset.title = x;
      }
    }
    const sample = new Sample();
    sample.componentMeta = {props: [
      {name: 'value', type: 'string'},
      {name: 'title', type: 'string'},
      {name: 'fixed', type: 'string'},
    ]};
    const props = new Map(sample.getProperties().map((p) => [p.name, p]));

    expect(props.get('value')!.get!(sample), 'a');
    props.get('value')!.set!(sample, 'b');
    expect(sample.value.value, 'b');

    expect(props.get('title')!.get!(sample), 'untitled');
    props.get('title')!.set!(sample, 'renamed');
    expect(sample.title, 'renamed');
    expect(sample.root.dataset.title, 'renamed');

    expect(props.get('fixed')!.get!(sample), 'construction-only');
    expect(props.get('fixed')!.set === undefined, true);
    sample.dispose();
  }, {owner: 'askalkin@datagrok.ai'});

  test('debug snapshot/diff', async () => {
    const debug = DG.U2.debug;
    const wasEnabled = debug.enabled;
    if (!wasEnabled)
      debug.enable({orphanSeconds: 0});
    try {
      const snapshot = debug.snapshot();
      const control = new DG.U2.Control();
      control.scope.listen(document, 'click', () => {});
      const report = debug.diff(snapshot);
      expect(report.scopes.length, 1);
      expect(report.listeners.length, 1);
      expect(report.total, 2);
      control.dispose();
      expect(debug.diff(snapshot).total, 0);
      expect(debug.liveCount() >= 0, true);
    } finally {
      if (!wasEnabled)
        debug.disable();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('DG.U2.state of a control', async () => {
    const control = new DG.U2.Control();
    control.root.dataset.u2 = 'sample';
    control.root.setAttribute('aria-label', 'Sample');
    const row = ui.divText('Alpha');
    row.dataset.u2Item = 'a';
    control.root.append(row);
    document.body.append(control.root);
    try {
      const state = DG.U2.state(row)!;
      expect(state.kind, 'sample');
      expect(state.label, 'Sample');
      expect(state.enabled, true);
      expect(state.visible, true);
      expect(DG.U2.status(control.root)!.state!.kind, 'sample');
      expect(await DG.U2.revealItem(control.root, 'alpha'), row);
      control.root.setAttribute('aria-disabled', 'true');
      control.visible = false;
      expect(DG.U2.state(control.root)!.enabled, false);
      expect(DG.U2.state(control.root)!.visible, false);
      expect(DG.U2.state(document.createElement('div')), null);
    } finally {
      control.root.remove();
      control.dispose();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('DG.Widget status carries state', async () => {
    const widget = DG.Widget.fromRoot(ui.divText('w'));
    widget.addStatusProvider('busy', () => ({state: {busy: true}}));
    try {
      const state = widget.getWidgetStatus().state!;
      expect(state.kind, widget.type);
      expect(state.enabled, true);
      expect(state.busy, true);
      expect(DG.U2.state(widget.root), null);
    } finally {
      widget.removeStatusProvider('busy');
      widget.detach();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('Widget.onEvent sees fireEvent', async () => {
    class Sample extends DG.Widget {
      raise(x: number): void {
        this.fireEvent('sample-raised', x);
      }
    }
    const widget = new Sample(ui.div());
    const seen: unknown[] = [];
    const sub = widget.onEvent('sample-raised').subscribe((x) => seen.push(x));
    const clicks = widget.onEvent('click').subscribe((e) => seen.push(e.type));
    widget.raise(1);
    widget.root.dispatchEvent(new MouseEvent('click'));
    sub.unsubscribe();
    clicks.unsubscribe();
    widget.raise(2);
    expect(seen.join(','), '1,click');
    widget.detach();
  }, {owner: 'askalkin@datagrok.ai'});

  test('Scope.bind after dispose', async () => {
    const scope = new DG.U2.Scope();
    let calls = 0;
    const handler = scope.bind(() => ++calls);
    expect(handler(), 1);
    scope.dispose();
    expect(handler() === undefined, true);
    expect(calls, 1);
  }, {owner: 'askalkin@datagrok.ai'});
});
