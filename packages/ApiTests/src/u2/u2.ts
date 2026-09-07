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
});
