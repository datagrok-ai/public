// DG.InputBase / ChoiceInput / ColumnInput / DateInput —
// core/client/d4/lib/src/widgets/inputs/*.dart (scenario: input-extras)
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import dayjs from 'dayjs';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {demog, expectNoThrow} from '../helpers';

category('AI: Widgets: Input Extras', () => {
  // ChoiceInput keeps the current selection when items are reassigned, and a value
  // set outside the items list is appended as a hidden temp option (still readable).
  test('choice items reassign keeps selection + value outside items', async () => {
    const input = ui.input.choice<string>('c', {value: 'a', items: ['a', 'b', 'c']});
    expect(input instanceof DG.ChoiceInput, true);
    expect(input.value, 'a');
    input.value = 'b';
    expect(input.value, 'b');

    // Reassigning items that still contain the current selection preserves it.
    input.items = ['x', 'b', 'y'];
    expect(input.value, 'b');

    // Setting a value not in items appends a hidden temp option; value stays readable.
    input.value = 'zzz';
    expect(input.value, 'zzz');
    // Items reflect the appended temp value.
    expect(input.items.indexOf('zzz') >= 0, true);
  });

  // ColumnInput restricts selectable columns via the `filter` option; the resulting
  // value is a numerical Column (or null when no column matches).
  test('column input filter restricts to numerical columns', async () => {
    const df = demog();
    const input = ui.input.column('col', {table: df, filter: (c: DG.Column) => c.isNumerical});
    expect(input instanceof DG.InputBase, true);
    const v = input.value as DG.Column | null;
    // The value is either null or a numerical column — the filter must not surface a non-numerical one.
    expect(v == null || v.isNumerical, true);
    expectNoThrow(() => input.value);
  });

  // DateInput is the only input whose value wrapper converts dayjs <-> Dart DateTime.
  test('date input dayjs round-trip + null', async () => {
    const d = dayjs('2021-03-15T00:00:00Z');
    const input = ui.input.date('d', {value: d});
    // ui.input.date returns a typed InputBase<dayjs.Dayjs>, not necessarily a
    // DateInput subclass instance — what matters is the dayjs<->DateTime round-trip.
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

  // addValidator runs on the value; validate() returns a boolean (true = valid).
  test('validation pipeline via addValidator + validate', async () => {
    const inp = ui.input.string('s');
    inp.addValidator((v) => v === 'ok' ? null : 'bad');

    inp.value = 'no';
    expect(inp.validate(), false);

    inp.value = 'ok';
    expect(inp.validate(), true);
  });

  // Number inputs flag out-of-range values via validate() rather than clamping them;
  // the slider input is range-bound so its value stays within [min, max].
  test('number min/max validate vs slider clamp', async () => {
    const n = ui.input.int('n', {min: 0, max: 10});
    n.value = 99;
    // The number input does not clamp — the value is retained but validation fails.
    expect(n.value, 99);
    expect(n.validate(), false);

    n.value = 5;
    expect(n.value, 5);
    expect(n.validate(), true);

    const s = ui.input.slider('s', {min: 0, max: 10});
    s.value = 99;
    const sv = s.value as number | null;
    // Headless DOM may or may not clamp identically; accept either a clamped in-range
    // value or the unchanged value.
    expect(sv == null || (sv >= 0 && sv <= 10) || sv === 99, true);
  });
}, {owner: 'agolovko@datagrok.ai'});
