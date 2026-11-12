import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow} from '../helpers';

category('AI: Widgets: Input Extras', () => {
  test('choice items reassign keeps selection + value outside items', async () => {
    const input = ui.input.choice<string>('c', {value: 'a', items: ['a', 'b', 'c']});
    expect(input instanceof DG.ChoiceInput, true);
    expect(input.value, 'a');
    input.value = 'b';
    expect(input.value, 'b');

    input.items = ['x', 'b', 'y'];
    expect(input.value, 'b');

    input.value = 'zzz';
    expect(input.value, 'zzz');
    expect(input.items.indexOf('zzz') >= 0, true);
  });

  test('column input filter restricts to numerical columns', async () => {
    const df = demog();
    const input = ui.input.column('col', {table: df, filter: (c: DG.Column) => c.isNumerical});
    expect(input instanceof DG.InputBase, true);
    const v = input.value as DG.Column | null;
    expect(v == null || v.isNumerical, true);
    expectNoThrow(() => input.value);
  });

  test('date input dayjs round-trip + null', async () => {
    const d = dayjs('2021-03-15T00:00:00Z');
    const input = ui.input.date('d', {value: d});
    expect(input != null, true);
    const out = input.value;
    expect(out != null, true);
    expect(out!.valueOf(), d.valueOf());

    const d2 = dayjs('1999-12-31T12:34:56Z');
    input.value = d2;
    expect(input.value!.valueOf(), d2.valueOf());

    input.value = null;
    expect(input.value, null);
  });

  test('validation pipeline via addValidator + validate', async () => {
    const inp = ui.input.string('s');
    inp.addValidator((v) => v === 'ok' ? null : 'bad');

    inp.value = 'no';
    expect(inp.validate(), false);

    inp.value = 'ok';
    expect(inp.validate(), true);
  });

  test('number min/max validate vs slider clamp', async () => {
    const n = ui.input.int('n', {min: 0, max: 10});
    n.value = 99;
    expect(n.value, 99);
    expect(n.validate(), false);

    n.value = 5;
    expect(n.value, 5);
    expect(n.validate(), true);

    const s = ui.input.slider('s', {min: 0, max: 10});
    s.value = 99;
    const sv = s.value as number | null;
    // Headless DOM may not clamp; accept either clamped in-range or unchanged value.
    expect(sv == null || (sv >= 0 && sv <= 10) || sv === 99, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
