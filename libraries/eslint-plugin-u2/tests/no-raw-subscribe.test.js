'use strict';

const tester = require('./tester');
const rule = require('../rules/no-raw-subscribe');

const U2 = `import {subscribe, div} from '@datagrok-libraries/u2';\n`;
const raw = (code) => ({code: U2 + code, errors: [{messageId: 'raw'}]});

tester.run('no-raw-subscribe', rule, {
  valid: [
    U2 + `subscribe(df.onCurrentRowChanged, () => refresh());`,
    U2 + `self.subscribe(grok.events.onViewRemoved, (v) => close(v));`,
    U2 + `host.subscribe(df.onValuesChanged, this.countMoved);`,
    U2 + `owner.subscribe(param.onChanged, edited);`,
    U2 + `ViewHost.of(view).subscribe(grok.events.onCurrentViewChanged, (v) => pick(v));`,
    U2 + `host.task(async (t) => { await load(); t.subscribe(df.onFilterChanged, refresh); });`,
    `this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => update()));`,
  ],
  invalid: [
    raw(`this.multiView.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
      if (grok.shell.v?.name === this.name) this.refresh();
    }));`),
    raw(`this._rowSub = df.onCurrentRowChanged.subscribe(() => {
      editor.load(df.currentRowIdx);
    });`),
    raw(`const changes = monomers.onChanged.subscribe((symbols) => picked.value = symbols.length);`),
    raw(`input.onChanged.subscribe(refresh);`),
    {code: U2 + `df.onValuesChanged.subscribe(() => recount());`, errors: [{message: 'A raw `.subscribe(…)` is ' +
      'owned by nothing and outlives the control unless unsubscribed by hand. Subscribe on the owner: ' +
      '`self.subscribe(source, handler)` in a component, `host.subscribe(source, handler)` on a view\'s `ViewHost`, ' +
      '`owner.subscribe(source, handler)` on any other owner. It ends with the owner and contains the handler.'}]},
    raw(`stream.subscribe({next: (x) => log(x), error: (e) => fail(e)});`),
    {code: `import * as DG from 'datagrok-api/dg';
      const c = new DG.U2.Control(div());
      grid.onCellClick.subscribe(() => c.dispose());`, errors: [{messageId: 'raw'}]},
  ],
});
