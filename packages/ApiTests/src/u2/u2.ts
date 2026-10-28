import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject} from 'rxjs';
import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';

category('U2', () => {
  let mode: DG.U2.UsageMode;

  before(async () => {
    mode = DG.U2.Usage.mode;
    DG.U2.Usage.mode = 'throw';
  });

  after(async () => {
    DG.U2.Usage.mode = mode;
  });

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
    const owner = DG.U2.Scope.root('api-tests debug');
    try {
      const snapshot = debug.snapshot();
      const control = owner.build(() => new DG.U2.Control());
      control.scope.listen(document, 'click', () => {});
      const report = debug.diff(snapshot);
      expect(report.scopes.length, 1);
      expect(report.listeners.length, 1);
      expect(report.total, 2);
      control.dispose();
      expect(debug.diff(snapshot).total, 0);
      expect(debug.liveCount() >= 0, true);
    } finally {
      owner.dispose();
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

  test('version', async () => {
    const [major, minor] = DG.U2.version.split('.').map(Number);
    expect(major > 1 || major === 1 && minor >= 2, true, `DG.U2.version ${DG.U2.version}`);
  }, {owner: 'askalkin@datagrok.ai'});

  test('Usage', async () => {
    const usage = DG.U2.Usage;
    const gone = DG.U2.Scope.root('api-tests session');
    gone.dispose();
    try {
      usage.clear();
      for (let i = 0; i < 2; i++) {
        let thrown: unknown;
        try {
          gone.effect(() => {});
        }
        catch (e) {
          thrown = e;
        }
        expect(thrown instanceof DG.U2.U2UsageError, true);
        expect((thrown as DG.U2.U2UsageError).code, 'owner-disposed');
      }
      const reports = usage.reports();
      expect(reports.length, 1);
      expect(reports[0].code, 'owner-disposed');
      expect(reports[0].count, 2, 'what was thrown is kept too, once per call site');
      expect(reports[0].message.includes('owner.task(async (t) =>'), true, reports[0].message);
      usage.mode = 'off';
      usage.clear();
      gone.listen(document, 'click', () => {});
      expect(usage.reports().length, 0);
    } finally {
      usage.mode = 'throw';
      usage.clear();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('owner verbs: root, build, task, subscribe', async () => {
    const live = DG.U2.Scope.liveCount;
    const root = DG.U2.Scope.root('api-tests verbs');
    const built = root.build(() => new DG.U2.Control());
    built.own(() => {});
    let late: DG.U2.Control | undefined;
    await root.task(async (t) => {
      await Promise.resolve();
      late = t.build(() => new DG.U2.Control());
    }).done;
    const started: number[] = [];
    let release!: () => void;
    const gate = new Promise<void>((resolve) => release = resolve);
    const runs = [1, 2, 3].map((n) => root.task(async () => {
      started.push(n);
      if (n === 1)
        await gate;
    }, {serial: 'save'}));
    release();
    await Promise.all(runs.map((r) => r.done));
    expect(started.join(','), '1,3');
    const seen: number[] = [];
    const stream = new Subject<number>();
    root.subscribe<number>(stream, (x) => seen.push(x));
    stream.next(1);
    let skipped: unknown = 'not reached';
    const pending = root.task(async (t) => {
      await Promise.resolve();
      skipped = t.build(() => 'built');
    });
    root.dispose();
    await pending.done;
    stream.next(2);
    expect(seen.join(','), '1');
    expect(stream.observers.length, 0, 'unsubscribed with its owner');
    expect(skipped === undefined, true);
    expect(built.scope.isDisposed, true);
    expect(late!.scope.isDisposed, true);
    expect(DG.U2.Scope.liveCount, live);
  }, {owner: 'askalkin@datagrok.ai'});

  test('lens and afterFlush', async () => {
    const person = DG.U2.signal({name: 'Ann', age: 41});
    const name = DG.U2.lens(person, (p) => p.name, (p, value) => ({...p, name: value}));
    expect(name instanceof DG.U2.Signal, true);
    expect(DG.U2.isWritableSignal(name), true);
    name.value = 'Bo';
    expect(person.value.name, 'Bo');
    expect(person.value.age, 41);
    const before = person.value;
    name.value = 'Bo';
    expect(person.value === before, true, 'the same value writes nothing');

    const rows = DG.U2.signal(1);
    const root = DG.U2.Scope.root('api-tests flush');
    let shown = 0;
    root.effect(() => {
      shown = rows.value;
    });
    const seen: number[] = [];
    DG.U2.batch(() => {
      rows.value = 2;
      DG.U2.afterFlush(() => seen.push(shown));
      seen.push(shown);
    });
    expect(seen.join(','), '1,2');
    root.dispose();
  }, {owner: 'askalkin@datagrok.ai'});

  test('operability: the six rules, on the old inert-rail rule and a real layout', async () => {
    const style = document.createElement('style');
    style.textContent = '.api-tests-u2-readonly .api-tests-u2-rail {pointer-events: none}';
    const row = (...children: HTMLElement[]) => {
      const el = ui.div(children);
      el.style.cssText = 'height: 120px; position: relative;';
      return el;
    };
    const iconButton = (label: string, title: string) => {
      const b = ui.button(ui.iconFA('pen'), () => {});
      b.setAttribute('aria-label', label);
      b.title = title;
      return b;
    };
    const rail = (...children: HTMLElement[]) => {
      const el = ui.div(children, 'api-tests-u2-rail');
      el.dataset.u2Part = 'options';
      el.style.cssText = 'height: 28px; display: flex;';
      return el;
    };
    const readOnly = ui.div([rail(iconButton('Edit Campaign Name', 'Edit Campaign Name'))], 'api-tests-u2-readonly');
    const hidden = ui.button('Under', () => {});
    const cover = ui.div([]);
    cover.style.cssText = 'position: absolute; left: 0; top: 0; width: 200px; height: 60px; background: white;';
    const pane = ui.div([ui.button('Run', () => {})]);
    pane.setAttribute('inert', '');
    const off = ui.button('Save', () => {});
    off.setAttribute('aria-disabled', 'true');
    const drawing = ui.div([]);
    drawing.style.cssText = 'width: 160px; height: 100px; flex-shrink: 0;';
    const root = ui.divV([row(readOnly), row(hidden, cover), row(pane), row(off), row(iconButton('Rename', 'Edit')),
      row(rail(drawing))]);
    root.style.cssText = 'position: fixed; left: 0; top: 0; width: 400px; background: white; z-index: 100000;';
    document.head.append(style);
    document.body.append(root);
    try {
      const found = DG.U2.operability(root, {hitTest: true});
      expect([...new Set(found.map((f) => f.rule))].sort().join(','), ['covered', 'disabled-without-reason', 'inert',
        'name-mismatch', 'no-pointer', 'oversized-option'].join(','), found.map((f) => f.detail).join('; '));
      expect(found.find((f) => f.rule === 'covered')?.element, hidden);
      expect(found.find((f) => f.rule === 'no-pointer')?.detail,
        'button "Edit Campaign Name" takes no pointer: pointer-events: none on div.api-tests-u2-rail');
      style.remove();
      off.setAttribute('aria-description', 'Name the campaign first');
      expect(DG.U2.operability(root).map((f) => f.rule).sort().join(','),
        'inert,name-mismatch,oversized-option', 'today\'s rail and a reason clear two; covered needs hitTest');
    } finally {
      style.remove();
      root.remove();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('component: a handle, a name, the owner verbs and disabledReason', async () => {
    const count = DG.U2.signal(0);
    const counter = DG.U2.component((self) => {
      const text = ui.divText('0');
      self.effect(() => text.textContent = String(count.value));
      return {root: text, get doubled() {
        return count.value * 2;
      }, increment: () => count.value++};
    }, {name: 'api-tests-counter'});
    counter.increment();
    expect(counter.doubled, 2);
    expect(counter.root.textContent, '1');
    expect(counter.root.dataset.u2Name, 'api-tests-counter');
    const owner = counter as unknown as DG.U2.Control & DG.U2.OwnerVerbs;
    const child = owner.build(() => new DG.U2.Control());
    child.own(() => {});
    counter.root.setAttribute('aria-disabled', 'true');
    counter.root.setAttribute('aria-description', 'Pick a campaign first');
    expect(counter.getState().disabledReason, 'Pick a campaign first');
    counter.dispose();
    count.value = 5;
    expect(counter.root.textContent, '1');
    expect(child.scope.isDisposed, true, 'what the control built goes with it');
    const alias = DG.U2.Control.build(() => ui.div(), {name: 'alias'});
    expect(alias.root.dataset.u2Name, 'alias', 'the deprecated Control.build, which published bundles call');
    alias.dispose();
  }, {owner: 'askalkin@datagrok.ai'});

  test('reserved-key: a library subclass is not reported; a field taking a control member throws', async () => {
    const debug = DG.U2.debug;
    const wasEnabled = debug.enabled;
    if (!wasEnabled)
      debug.enable({orphanSeconds: 0});
    class LibraryCard extends DG.Widget {}
    class ConsumerCard extends LibraryCard {}
    class EffectCard extends DG.Widget {
      constructor() {
        super(ui.div());
        Object.defineProperty(this, 'effect', {value: 1, enumerable: true});
      }
    }
    // the fields are checked a microtask after the constructor, so the throw reaches the page as an uncaught error
    const thrown: DG.U2.U2UsageError[] = [];
    const caught = (e: ErrorEvent) => {
      if (e.error instanceof DG.U2.U2UsageError) {
        thrown.push(e.error);
        e.preventDefault();
        e.stopImmediatePropagation();
      }
    };
    window.addEventListener('error', caught, true);
    try {
      new ConsumerCard(ui.div()).detach();
      await Promise.resolve();
      expect(thrown.length, 0);
      new EffectCard().detach();
      await Promise.resolve();
      expect(thrown.map((e) => e.code).join(','), 'reserved-key');
      // a release build minifies the class name: the member is what the message must name
      expect(thrown[0].message.includes('.effect overrides Control\'s own `effect()`'), true, thrown[0].message);
    } finally {
      window.removeEventListener('error', caught, true);
      grok.shell.clearLastError();
      DG.U2.Usage.clear();
      if (!wasEnabled)
        debug.disable();
    }
  }, {owner: 'askalkin@datagrok.ai'});

  test('Scope.thenable', async () => {
    const promise = Promise.resolve(1);
    const custom = {then: () => {}};
    const row = DG.DataFrame.fromColumns([DG.Column.fromStrings('name', ['a'])]).row(0);
    const throwing = new Proxy({}, {get: () => {
      throw new Error('No such column');
    }});
    expect(DG.U2.Scope.thenable(promise), promise);
    expect(DG.U2.Scope.thenable(custom), custom);
    expect(DG.U2.Scope.thenable(row) === undefined, true);
    expect(DG.U2.Scope.thenable(throwing) === undefined, true);
    expect(DG.U2.Scope.thenable('then') === undefined, true);

    const scope = new DG.U2.Scope();
    const errors: unknown[] = [];
    scope.onError = (e) => errors.push(e);
    expect(scope.bind(() => row)(), row);
    expect(scope.bind(() => throwing)() === throwing, true);
    expect(errors.length, 0);
    scope.dispose();
  }, {owner: 'askalkin@datagrok.ai'});

  test('a static control holds no live scope', async () => {
    const owner = DG.U2.Scope.root('api-tests static');
    const base = DG.U2.Scope.liveCount;
    const control = owner.build(() => new DG.U2.Control());
    control.scope.run(() => control.root.append(ui.divText('static')));
    control.name = 'static-control';
    expect(DG.U2.Scope.liveCount, base);
    control.own(() => {});
    expect(DG.U2.Scope.liveCount, base + 1);
    control.dispose();
    expect(DG.U2.Scope.liveCount, base);
    owner.dispose();
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
