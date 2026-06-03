import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, after, awaitCheck} from '@datagrok-libraries/test/src/test';
import {closeView, awaitWebComponents} from './utils';

/** SimpleInputs2 has params: int a, double b, string c — the second input is the float. */
const FLOAT_INPUT_INDEX = 1;

async function openSimpleInputs(b: number): Promise<{view: DG.ViewBase, call: DG.FuncCall}> {
  await awaitWebComponents();
  const func = DG.Func.byName('Compute2:SimpleInputs2');
  const call = func.prepare({a: 1, b, c: 'x'});
  const view = await grok.functions.call(
    'Compute2:RichFunctionViewEditor', {call},
  ) as unknown as DG.ViewBase;
  return {view, call};
}

async function getFloatInput(view: DG.ViewBase): Promise<HTMLInputElement> {
  await awaitCheck(
    () => view.root.querySelectorAll('dg-input-form input').length >= 3,
    'Form inputs not rendered', 15000,
  );
  const inputs = Array.from(view.root.querySelectorAll('dg-input-form input')) as HTMLInputElement[];
  return inputs[FLOAT_INPUT_INDEX];
}

function simulateTyping(input: HTMLInputElement, value: string): void {
  input.focus();
  const setter = Object.getOwnPropertyDescriptor(HTMLInputElement.prototype, 'value')!.set!;
  setter.call(input, value);
  input.dispatchEvent(new Event('input', {bubbles: true}));
  input.dispatchEvent(new Event('change', {bubbles: true}));
  input.blur();
}

category('InputForm: Float32 display cleanup', () => {
  let view: DG.ViewBase | undefined;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Cleans Float32 noise on bind: 0.2', async () => {
    const noisy = Math.fround(0.2); // 0.20000000298023224
    const result = await openSimpleInputs(noisy);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '0.2',
      `Expected "0.2", got "${bInput.value}"`, 5000,
    );
  });

  test('Cleans Float32 noise on bind: negative -0.2', async () => {
    const noisy = Math.fround(-0.2);
    const result = await openSimpleInputs(noisy);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '-0.2',
      `Expected "-0.2", got "${bInput.value}"`, 5000,
    );
  });

  test('Cleans Float32 noise on bind: 0.1', async () => {
    const noisy = Math.fround(0.1);
    const result = await openSimpleInputs(noisy);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '0.1',
      `Expected "0.1", got "${bInput.value}"`, 5000,
    );
  });

  test('Idempotent for clean short value: 0.5', async () => {
    const result = await openSimpleInputs(0.5);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '0.5',
      `Expected "0.5", got "${bInput.value}"`, 5000,
    );
  });

  test('Idempotent for integer-like value: 1', async () => {
    const result = await openSimpleInputs(1);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '1',
      `Expected "1", got "${bInput.value}"`, 5000,
    );
  });

  test('Handles modulus > 1: 1234.45345', async () => {
    const result = await openSimpleInputs(1234.45345);
    view = result.view;
    const bInput = await getFloatInput(view);
    // Float32(1234.45345) ≈ 1234.45349... — display via #0.### should be a
    // short form that parses back close to original.
    await awaitCheck(
      () => {
        const parsed = parseFloat(bInput.value);
        return !isNaN(parsed) && Math.abs(parsed - 1234.45345) < 0.01;
      },
      `Expected display close to 1234.45345, got "${bInput.value}"`, 5000,
    );
    // Should not contain trailing-noise digits like "...000000".
    expect(bInput.value.length <= 12, true,
      `Display should be short, got "${bInput.value}"`);
  });

  test('Handles large magnitude: 10000', async () => {
    const result = await openSimpleInputs(10000);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => parseFloat(bInput.value) === 10000,
      `Expected 10000, got "${bInput.value}"`, 5000,
    );
  });

  test('Handles small magnitude: 0.0001', async () => {
    const result = await openSimpleInputs(0.0001);
    view = result.view;
    const bInput = await getFloatInput(view);
    // The default mask #0.### keeps 3 decimals, so 0.0001 is below display
    // resolution and cleanly rounds to "0" — no scientific notation, no noise.
    await awaitCheck(
      () => bInput.value === '0',
      `Expected "0" (0.0001 rounds away at 3 decimals), got "${bInput.value}"`, 5000,
    );
  });

  test('Display stays clean after typing a noisy value', async () => {
    const result = await openSimpleInputs(0);
    view = result.view;
    const bInput = await getFloatInput(view);
    await awaitCheck(
      () => bInput.value === '0',
      'Initial value should be 0', 5000,
    );
    // Type a Float32-noisy string. Format on input keeps display clean after parse.
    simulateTyping(bInput, '0.20000000298023224');
    await awaitCheck(
      () => bInput.value === '0.2' || bInput.value === '0.20000000298023224',
      `After typing noisy value, expected clean display or as-typed, got "${bInput.value}"`,
      5000,
    );
  });

  test('Two consecutive forms both display cleanly', async () => {
    const first = await openSimpleInputs(Math.fround(0.2));
    view = first.view;
    const firstInput = await getFloatInput(view);
    await awaitCheck(
      () => firstInput.value === '0.2',
      `First form: expected "0.2", got "${firstInput.value}"`, 5000,
    );
    closeView(view);
    view = undefined;

    const second = await openSimpleInputs(Math.fround(0.3));
    view = second.view;
    const secondInput = await getFloatInput(view);
    await awaitCheck(
      () => secondInput.value === '0.3',
      `Second form: expected "0.3", got "${secondInput.value}"`, 5000,
    );
  });

  // Last — mutates the cached Property.format, which would bleed into earlier tests.
  test('Pre-existing input format is preserved (not overridden)', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:SimpleInputs2');
    const call = func.prepare({a: 1, b: 3.14159, c: 'x'});
    const bParam = call.inputParams['b'];
    expect(bParam != null, true, 'b inputParam should exist');
    const originalFormat = (bParam.property as any).format;
    (bParam.property as any).format = '#0.00';
    try {
      view = await grok.functions.call(
        'Compute2:RichFunctionViewEditor', {call},
      ) as unknown as DG.ViewBase;
      const bInput = await getFloatInput(view);
      await awaitCheck(
        () => bInput.value === '3.14',
        `Pre-existing format should win, expected "3.14", got "${bInput.value}"`, 5000,
      );
    } finally {
      (bParam.property as any).format = originalFormat ?? '';
    }
  });
});
