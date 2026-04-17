import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, test, expect, after, awaitCheck} from '@datagrok-libraries/test/src/test';
import {closeView, awaitWebComponents} from './utils';

/** Find and click a button by exact label text. */
async function clickBigButton(root: HTMLElement, label: string, timeout = 10000): Promise<void> {
  await awaitCheck(() => {
    return Array.from(root.querySelectorAll('button')).some(
      (b) => b.textContent?.trim() === label,
    );
  }, `Button "${label}" not found`, timeout);
  const btn = Array.from(root.querySelectorAll('button')).find(
    (b) => b.textContent?.trim() === label,
  )!;
  btn.click();
}

/** Wait for Run button to disappear (signals run completion in RFV). */
async function awaitRunComplete(root: HTMLElement, timeout = 30000): Promise<void> {
  await awaitCheck(() => {
    return !Array.from(root.querySelectorAll('button')).some(
      (b) => b.textContent?.trim() === 'Run',
    );
  }, 'Run did not complete (Run button still visible)', timeout);
}

// ---- RFV Editor tests ----

category('Editors: RFV smoke', () => {
  let view: DG.ViewBase | undefined;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Opens for SimpleInputs2 with form and Run button', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:SimpleInputs2');
    const call = func.prepare();
    view = await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelector('dg-input-form') !== null,
      'RFV should contain dg-input-form', 15000,
    );
    await awaitCheck(
      () => view!.root.querySelectorAll('dg-input-form input').length >= 3,
      'SimpleInputs2 has 3 params — form should have at least 3 inputs', 5000,
    );
    await awaitCheck(() => {
      return Array.from(view!.root.querySelectorAll('button')).some(
        (b) => b.textContent?.trim() === 'Run',
      );
    }, 'Run button should be visible initially', 5000);
    await awaitCheck(
      () => view!.root.querySelector('[dock-spawn-title="Inputs"]') !== null,
      'Inputs dock panel not rendered', 5000,
    );
  });

  test('Opens for ObjectCooling2 with grouped inputs', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:ObjectCooling2');
    const call = func.prepare();
    view = await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelectorAll('dg-input-form input').length >= 5,
      'ObjectCooling2 form should have at least 5 input fields', 15000,
    );
    await awaitCheck(() => {
      const text = view!.root.textContent ?? '';
      return text.includes('Ambient') || text.includes('Initial') || text.includes('Surface');
    }, 'ObjectCooling2 form should show input labels', 5000);
  });
});

// ---- RFV run tests ----

category('RFV: run SimpleInputs2', () => {
  let view: DG.ViewBase | undefined;
  let call: DG.FuncCall;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Fill inputs, run, check outputs in DOM', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:SimpleInputs2');
    call = func.prepare({a: 42, b: 3.14, c: 'test-value'});
    view = await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelectorAll('input').length >= 3,
      'Form inputs not rendered', 15000,
    );

    await clickBigButton(view.root, 'Run');
    await awaitRunComplete(view.root);

    expect(call.outputs.get('aout'), 42, 'aout should be 42');
    expect(call.outputs.get('bout'), 3.14, 'bout should be 3.14');
    expect(call.outputs.get('cout'), 'test-value', 'cout should be test-value');

    await awaitCheck(() => {
      const text = view!.root.textContent ?? '';
      return text.includes('42') && text.includes('3.14') && text.includes('test-value');
    }, 'Scalar outputs not rendered in DOM', 10000);
  });
});

category('RFV: run ObjectCooling2', () => {
  let view: DG.ViewBase | undefined;
  let call: DG.FuncCall;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Fill inputs, run, verify outputs', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:ObjectCooling2');
    call = func.prepare({
      ambTemp: 20, initTemp: 80, desiredTemp: 30,
      area: 0.06, heatCap: 4200, heatTransferCoeff: 8.3, simTime: 100,
    });
    view = await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelectorAll('input').length >= 3,
      'Form inputs not rendered', 15000,
    );

    await clickBigButton(view.root, 'Run');
    await awaitRunComplete(view.root);

    const tempDiff = call.outputs.get('tempDiff');
    expect(tempDiff, 60, 'tempDiff should be initTemp - ambTemp = 60');

    const coolingFactor = call.outputs.get('coolingFactor');
    expect(typeof coolingFactor === 'number' && coolingFactor > 0, true,
      `coolingFactor should be positive, got ${coolingFactor}`);

    const sim = call.outputs.get('simulation') as DG.DataFrame;
    expect(sim instanceof DG.DataFrame, true, 'simulation should be a DataFrame');
    const tempCol = sim.col('Temperature');
    expect(tempCol !== null, true, 'simulation should have Temperature column');
    if (tempCol && tempCol.length > 0) {
      const firstTemp = tempCol.get(0);
      expect(Math.abs(firstTemp - 80) < 1, true,
        `First temperature should be ~80, got ${firstTemp}`);
    }

    await awaitCheck(() => {
      const text = view!.root.textContent ?? '';
      return text.includes('60');
    }, 'Scalar output (tempDiff=60) not rendered in DOM', 10000);
  });
});

// ---- Pipeline tests ----

category('Pipeline: MockPipeline1 static', () => {
  let view: DG.ViewBase | undefined;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Renders tree with steps and Proceed card', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:MockPipeline1');
    const call = func.prepare();
    view = await grok.functions.call(
      'Compute2:TreeWizardEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelector('.mtl-tree') !== null,
      'Tree not rendered', 15000,
    );

    // Verify tree nodes
    await awaitCheck(() => {
      const labels = view!.root.querySelectorAll('.mtl-tree .mtl-ml');
      const texts = Array.from(labels).map((l) => l.textContent?.trim());
      return texts.includes('step1') && texts.includes('step2');
    }, 'Tree should show step1 and step2', 10000);

    // Verify pipeline view shows Proceed card
    await awaitCheck(() => {
      return view!.root.textContent?.includes('Proceed to the first step') ?? false;
    }, 'Should show Proceed to the first step card', 10000);
  });
});

category('Pipeline: MockPipeline2 sequential', () => {
  let view: DG.ViewBase | undefined;

  after(async () => {
    if (view) closeView(view);
    view = undefined;
  });

  test('Renders tree with initial steps and Proceed card', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:MockPipeline2');
    const call = func.prepare();
    view = await grok.functions.call(
      'Compute2:TreeWizardEditor', {call},
    ) as unknown as DG.ViewBase;

    await awaitCheck(
      () => view!.root.querySelector('.mtl-tree') !== null,
      'Tree not rendered', 15000,
    );

    // Verify tree shows initial steps (add, mul, cooling)
    await awaitCheck(() => {
      const labels = view!.root.querySelectorAll('.mtl-tree .mtl-ml');
      const texts = Array.from(labels).map((l) => l.textContent?.trim());
      const count = (texts.some((t) => t === 'add') ? 1 : 0)
        + (texts.some((t) => t === 'mul') ? 1 : 0)
        + (texts.some((t) => t?.includes('cooling')) ? 1 : 0);
      return count >= 2;
    }, 'Tree should show at least 2 of: add, mul, cooling', 10000);

    // Verify pipeline view shows Proceed card
    await awaitCheck(() => {
      return view!.root.textContent?.includes('Proceed to the first step') ?? false;
    }, 'Should show Proceed to the first step card', 10000);
  });
});
