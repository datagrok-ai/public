import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as Vue from 'vue';
import {category, test, expect, waitUntil, awaitCheck, delay} from '@datagrok-libraries/test/src/test';
import {OUTPUT_OUTDATED_PATH} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {historyUtils} from '@datagrok-libraries/compute-utils/history-utils/src/history-utils';
import {closeViewsAfter, awaitWebComponents, deleteSavedWorkflow} from './utils';
import {DAY_CALCULATION, uniqueName, withRunStorage} from './run-storage';

/** Find and click a button by exact label text. */
async function clickBigButton(root: HTMLElement, label: string, timeout = 10000): Promise<void> {
  await waitUntil(() => {
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
  await waitUntil(() => {
    return !Array.from(root.querySelectorAll('button')).some(
      (b) => b.textContent?.trim() === 'Run',
    );
  }, 'Run did not complete (Run button still visible)', timeout);
}

// ---- RFV Editor tests ----

category('Editors: RFV smoke', () => {
  const track = closeViewsAfter();

  test('Opens for SimpleInputs2 with form and Run button', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:SimpleInputs2');
    const call = func.prepare();
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelector('dg-input-form') !== null, 'RFV should contain dg-input-form', 15000);
    await waitUntil(() => view.root.querySelectorAll('dg-input-form input').length >= 3,
      'SimpleInputs2 has 3 params — form should have at least 3 inputs', 5000);
    await waitUntil(() => {
      return Array.from(view.root.querySelectorAll('button')).some(
        (b) => b.textContent?.trim() === 'Run',
      );
    }, 'Run button should be visible initially', 5000);
    await waitUntil(() => view.root.querySelector('[dock-spawn-title="Inputs"]') !== null,
      'Inputs dock panel not rendered', 5000);
  });

  test('Opens for ObjectCooling2 with grouped inputs', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:ObjectCooling2');
    const call = func.prepare();
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelectorAll('dg-input-form input').length >= 5,
      'ObjectCooling2 form should have at least 5 input fields', 15000);
    await waitUntil(() => {
      const text = view.root.textContent ?? '';
      return text.includes('Ambient') || text.includes('Initial') || text.includes('Surface');
    }, 'ObjectCooling2 form should show input labels', 5000);
  });
});

// ---- RFV run tests ----

category('RFV: run SimpleInputs2', () => {
  const track = closeViewsAfter();
  let call: DG.FuncCall;


  test('Fill inputs, run, check outputs in DOM', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:SimpleInputs2');
    call = func.prepare({a: 42, b: 3.14, c: 'test-value'});
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelectorAll('input').length >= 3, 'Form inputs not rendered', 15000);
    await waitUntil(() => ribbonIcons(view, 'exchange').length === 1, 'Expected one compare icon', 5000);

    await clickBigButton(view.root, 'Run');
    await awaitRunComplete(view.root);

    expect(call.outputs.get('aout'), 42, 'aout should be 42');
    expect(call.outputs.get('bout'), 3.14, 'bout should be 3.14');
    expect(call.outputs.get('cout'), 'test-value', 'cout should be test-value');

    await waitUntil(() => {
      const text = view.root.textContent ?? '';
      return text.includes('42') && text.includes('3.14') && text.includes('test-value');
    }, 'Scalar outputs not rendered in DOM', 10000);
  });
});

category('RFV: viewers hook', () => {
  const track = closeViewsAfter();

  test('Opens with the viewers hook made', async () => {
    await awaitWebComponents();
    let made = false;
    const applied: {ioName: string, type: string}[] = [];
    const subs = [
      grok.events.onCustomEvent('compute2-test-viewers-hook-made').subscribe(() => made = true),
      grok.events.onCustomEvent('compute2-test-viewers-hook-applied').subscribe((e) => applied.push(e)),
    ];
    try {
      const call = DG.Func.byName('Compute2:ViewersHookModel').prepare();
      const view = track(await grok.functions.call(
        'Compute2:RichFunctionViewEditor', {call},
      ) as unknown as DG.ViewBase);
      expect(made, true, 'viewers hook should be made before the view opens');

      await clickBigButton(view.root, 'Run');
      await awaitRunComplete(view.root);
      await waitUntil(() => applied.some((e) => e.ioName === 'df' && e.type === 'Grid'),
        'viewers hook not applied to the df grid', 15000);
    } finally {
      subs.forEach((sub) => sub.unsubscribe());
    }
  });
});

category('RFV: run ObjectCooling2', () => {
  const track = closeViewsAfter();
  let call: DG.FuncCall;


  test('Fill inputs, run, verify outputs', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:ObjectCooling2');
    call = func.prepare({
      ambTemp: 20, initTemp: 80, desiredTemp: 30,
      area: 0.06, heatCap: 4200, heatTransferCoeff: 8.3, simTime: 100,
    });
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelectorAll('input').length >= 3, 'Form inputs not rendered', 15000);

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

    await waitUntil(() => {
      const text = view.root.textContent ?? '';
      return text.includes('60');
    }, 'Scalar output (tempDiff=60) not rendered in DOM', 10000);
  });
});

category('RFV: core validation', () => {
  const track = closeViewsAfter();

  const open = async (preset: Record<string, any>) => {
    await awaitWebComponents();
    const call = DG.Script.create(`//name: RfvCoreValidation
//language: javascript
//input: int ranged = 5 {min: 0; max: 10}
//input: string code = "12ab" {validator: startsWith(value, "12")}
//input: string req
//output: int res
res = ranged;`).prepare(preset);
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);
    // the form reports its validity to the view before the instance becomes readable here
    const formOf = () => (view.root.querySelector('dg-input-form') as any)?.formInst as DG.InputForm | undefined;
    await waitUntil(() => !!formOf(), 'Form not rendered', 15000);
    const runDisabled = () => Array.from(view.root.querySelectorAll('button'))
      .find((b) => b.textContent?.trim() === 'Run')!.getAttribute('aria-disabled') === 'true';
    return {view, call, form: formOf()!, runDisabled};
  };

  test('Opened with valid values can run', async () => {
    const {form, runDisabled} = await open({ranged: 5, code: '12ab', req: 'x'});
    expect(form.isValid, true);
    expect(runDisabled(), false);
  });

  test('Opened with invalid values blocks the run until fixed', async () => {
    const {view, call, form, runDisabled} = await open({ranged: 50, code: 'xyz', req: 'x'});
    expect(form.isValid, false);
    expect(runDisabled(), true);

    form.getInput('ranged').value = 3;
    form.getInput('code').value = '12cd';
    await waitUntil(() => !runDisabled(), 'Run should be enabled once inputs are fixed', 5000);
    await clickBigButton(view.root, 'Run');
    await awaitRunComplete(view.root);
    expect(call.outputs.get('res'), 3);
  });

  test('Opened with an empty required input blocks the run', async () => {
    const {form, runDisabled} = await open({ranged: 5, code: '12ab'});
    expect(form.isValid, false);
    expect(runDisabled(), true);
  });
});

category('RFV: core form layout', () => {
  const track = closeViewsAfter();

  const open = async () => {
    await awaitWebComponents();
    const call = DG.Script.create(`//name: RfvCoreFormLayout
//language: javascript
//input: bool advanced = false
//input: double rate = 1.5 {visible: advanced}
//input: int steps = 10 {visible: advanced}
//input: double basic = 3 {visible: !advanced}
//input: double always = 7
//input: string pick = "a" {choices: ["a", "b", "c"]}
//input: list<string> tags = ["alpha"] {choices: ["alpha", "beta", "gamma"]}
//output: int res
res = steps;`).prepare();
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);
    grok.shell.addView(view);
    const formOf = () => (view.root.querySelector('dg-input-form') as any)?.formInst as DG.InputForm | undefined;
    await waitUntil(() => !!formOf(), 'Form not rendered', 15000);
    return formOf()!;
  };

  test('The view keeps inputs hidden by a visible expression on open', async () => {
    // the core hides inputs bypassing this setter, so only the view's writes land here
    const shown: string[] = [];
    const proto = DG.InputBase.prototype;
    const visible = Object.getOwnPropertyDescriptor(proto, 'visible')!;
    Object.defineProperty(proto, 'visible', {...visible, set(this: DG.InputBase, value: boolean) {
      if (value)
        shown.push(this.property?.name);
      visible.set!.call(this, value);
    }});
    try {
      const form = await open();
      // status icons are injected in the same watcher flush as the view's visibility
      await waitUntil(() => !!form.getInput('rate').root.querySelector('.rfv2-validation-icon'),
        'Status icons not injected', 5000);
    } finally {
      Object.defineProperty(proto, 'visible', visible);
    }
    expect(shown.join(), '', 'The view showed inputs the core hid');
  });

  test('A multiple choice input in a condensed form shows its list', async () => {
    const form = await open();
    await waitUntil(() => form.root.classList.contains('ui-form-condensed'), 'The form is not condensed', 5000);
    const list = form.getInput('tags').root.querySelector('.ui-input-editor')!;
    await waitUntil(() => list.getBoundingClientRect().width > 0, 'The list has no width', 5000);
  });

  test('Choice inputs in a condensed form leave no more space below than other inputs', async () => {
    const form = await open();
    await waitUntil(() => form.root.classList.contains('ui-form-condensed'), 'The form is not condensed', 5000);
    const spaceBelow = (name: string) => form.getInput(name).root.getBoundingClientRect().bottom -
      form.getInput(name).root.querySelector('.ui-input-editor')!.getBoundingClientRect().bottom;
    await waitUntil(() => spaceBelow('always') > 0, 'Inputs not laid out', 5000);
    expect(spaceBelow('pick'), spaceBelow('always'), 'Space below the select');
    expect(spaceBelow('tags'), spaceBelow('always'), 'Space below the multiple choice list');
  });
});

category('RFV: autorun', () => {
  const track = closeViewsAfter();

  /** Opens an autorun script, then sets `a` to 10 and, while the run that starts is in progress, to [next]. */
  const editDuringRun = async (body: string, next: number) => {
    await awaitWebComponents();
    const call = DG.Script.create(`//name: RfvAutorun
//language: javascript
//input: double a = 1 {min: 0}
//input: double b = 2
//output: double c
//meta.runOnInput: true
await new Promise((r) => setTimeout(r, 1500));
${body}`).prepare({a: 1, b: 2});
    const view = track(await grok.functions.call(
      'Compute2:RichFunctionViewEditor', {call},
    ) as unknown as DG.ViewBase);
    const formOf = () => (view.root.querySelector('dg-input-form') as any)?.formInst as DG.InputForm | undefined;
    await waitUntil(() => !!formOf(), 'Form not rendered', 15000);
    const isRunning = () => !!view.root.querySelector('.c2-compositor-overlay--active');
    await waitUntil(() => !isRunning(), 'Initial autorun did not finish', 10000);

    formOf()!.getInput('a').value = 10;
    await waitUntil(isRunning, 'Autorun did not start', 5000);
    formOf()!.getInput('a').value = next;
    return {call, isRunning};
  };

  test('An input changed during a run is run again', async () => {
    const {call} = await editDuringRun('c = a + b;', 40);
    await waitUntil(() => call.outputs.get('c') === 42, () => `Expected c = 42, got ${call.outputs.get('c')}`, 10000);
  });

  test('An input changed during a failed run is run again', async () => {
    const {call} = await editDuringRun('if (a === 10) throw new Error("Expected failure");\nc = a + b;', 40);
    await waitUntil(() => call.outputs.get('c') === 42, () => `Expected c = 42, got ${call.outputs.get('c')}`, 10000);
  });

  test('An input made invalid during a run leaves the outputs outdated and is not run', async () => {
    const {call, isRunning} = await editDuringRun('c = a + b;', -5);
    await waitUntil(() => !isRunning(), 'Run did not finish', 10000);
    await delay(1000);
    expect(isRunning(), false, 'An invalid input should not start a run');
    expect(call.outputs.get('c'), 12);
    expect(call.options[OUTPUT_OUTDATED_PATH], 'true', 'Outputs of a=10 should be outdated for a=-5');
  });
});

// ---- Pipeline tests ----

category('Pipeline: MockPipeline1 static', () => {
  const track = closeViewsAfter();

  test('A never saved workflow shows its root title until it is saved under its own name', async () => {
    await awaitWebComponents();
    const call = DG.Func.byName('Compute2:MockPipeline1').prepare();
    const view = track(await grok.functions.call(
      'Compute2:TreeWizardEditor', {call},
    ) as unknown as DG.ViewBase);
    const model = call.func.friendlyName;
    await waitUntil(() => view.name === `${model} - Pipeline 1 title`, 'View should show the root title', 5000);

    await waitUntil(() => ribbonIcons(view, 'save').length > 0, 'Expected a save icon', 5000);
    (ribbonIcons(view, 'save')[0] as HTMLElement).click();
    await waitUntil(() => document.querySelector('.d4-dialog .c2-runs-field input') !== null,
      'Save dialog not shown', 5000);
    const name = uniqueName('Own name');
    const input = document.querySelector<HTMLInputElement>('.d4-dialog .c2-runs-field input')!;
    expect(input.value, 'Pipeline 1 title');
    input.value = name;
    input.dispatchEvent(new Event('input'));
    await Vue.nextTick();
    Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog button'))
      .find((b) => b.textContent?.trim() === 'Save')!.click();
    await waitUntil(() => savedRunId(view) != null, 'URL should carry the saved run id', 5000);
    try {
      await waitUntil(() => view.name === `${model} - ${name}`, 'View should show the saved name', 5000);
    } finally {
      await deleteSavedWorkflow(savedRunId(view)!);
    }
  });

  test('Renders tree with steps and Next button', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:MockPipeline1');
    const call = func.prepare();
    const view = track(await grok.functions.call(
      'Compute2:TreeWizardEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelector('.mtl-tree') !== null, 'Tree not rendered', 15000);

    // Verify tree nodes
    await waitUntil(() => {
      const labels = view.root.querySelectorAll('.mtl-tree .mtl-ml');
      const texts = Array.from(labels).map((l) => l.textContent?.trim());
      return texts.includes('step1') && texts.includes('step2');
    }, 'Tree should show step1 and step2', 10000);
    await waitUntil(() => ribbonIcons(view, 'exchange').length === 1, 'Expected one compare icon', 5000);

    // Pipeline view's bottom bar shows the unified Next button (same as RFV).
    await waitUntil(() => {
      return Array.from(view.root.querySelectorAll('button')).some(
        (b) => b.textContent?.trim() === 'Next',
      );
    }, 'Should show Next button', 10000);
  });
});

category('Pipeline: MockPipeline2 sequential', () => {
  const track = closeViewsAfter();

  test('Renders tree with initial steps and Next button', async () => {
    await awaitWebComponents();
    const func = DG.Func.byName('Compute2:MockPipeline2');
    const call = func.prepare();
    const view = track(await grok.functions.call(
      'Compute2:TreeWizardEditor', {call},
    ) as unknown as DG.ViewBase);

    await waitUntil(() => view.root.querySelector('.mtl-tree') !== null, 'Tree not rendered', 15000);

    // Verify tree shows initial steps (add, mul, cooling)
    await waitUntil(() => {
      const labels = view.root.querySelectorAll('.mtl-tree .mtl-ml');
      const texts = Array.from(labels).map((l) => l.textContent?.trim());
      const count = (texts.some((t) => t === 'add') ? 1 : 0) +
        (texts.some((t) => t === 'mul') ? 1 : 0) +
        (texts.some((t) => t?.includes('cooling')) ? 1 : 0);
      return count >= 2;
    }, 'Tree should show at least 2 of: add, mul, cooling', 10000);

    // Pipeline view's bottom bar shows the unified Next button (same as RFV).
    await waitUntil(() => {
      return Array.from(view.root.querySelectorAll('button')).some(
        (b) => b.textContent?.trim() === 'Next',
      );
    }, 'Should show Next button', 10000);
  });
});

// ---- Compact single-step workflow ----

const hasButton = (root: HTMLElement, label: string) =>
  Array.from(root.querySelectorAll('button')).some((b) => b.textContent?.trim() === label);

// the editor view is not added to the shell here, so read the ribbon through the API
export const ribbonIcons = (view: DG.ViewBase, name: string) => view.getRibbonPanels().flat()
  .flatMap((el) => [el, ...Array.from(el.querySelectorAll('*'))])
  .filter((el) => el.classList.contains(`fa-${name}`));

export const savedRunId = (view: DG.ViewBase) => new URLSearchParams((view.path ?? '').split('?')[1] ?? '').get('id');

async function openCompact(nqName: string): Promise<DG.ViewBase> {
  await awaitWebComponents();
  const call = DG.Func.byName(nqName).prepare();
  const view = await grok.functions.call('Compute2:TreeWizardEditor', {call}) as unknown as DG.ViewBase;
  await waitUntil(() => view.root.querySelector('dg-input-form') !== null, 'Step form not rendered', 20000);
  await waitUntil(() => hasButton(view.root, 'Run'), 'Run button not rendered', 10000);
  return view;
}

for (const nqName of ['Compute2:MockSingleStepPipeline', 'Compute2:MockSingleStepNested']) {
  category(`Pipeline: single step compact (${nqName.split(':')[1]})`, () => {
    const track = closeViewsAfter();

    test('No tree, no navigation, one save and one share control', async () => {
      const view = track(await openCompact(nqName));
      expect(view.root.querySelector('.mtl-tree'), null, 'Steps tree should be hidden');
      expect(hasButton(view.root, 'Back'), false, 'Back should be hidden');
      expect(hasButton(view.root, 'Next'), false, 'Next should be hidden');
      expect(hasButton(view.root, 'Double initial temperature'), true, 'Step button action missing');
      expect(view.root.querySelector('[dock-spawn-title="Inputs"]') !== null, true, 'Inputs panel missing');
      expect(ribbonIcons(view, 'save').length, 1, 'Expected exactly one save icon');
      expect(ribbonIcons(view, 'exchange').length, 1, 'Expected exactly one compare icon');
      expect(ribbonIcons(view, 'folder-tree').length, 0, 'Tree toggle should be hidden');
      expect(ribbonIcons(view, 'share-alt').length <= 1, true, 'Expected at most one share icon');
    });

    test('Run, then save with id in the URL', async () => {
      const view = track(await openCompact(nqName));
      expect(ribbonIcons(view, 'save').length, 1, 'Expected exactly one save icon before run');

      await clickBigButton(view.root, 'Run');
      await awaitRunComplete(view.root);
      await waitUntil(() => (view.root.textContent ?? '').includes('80'), 'tempDiff (100 - 20) not rendered', 10000);

      (ribbonIcons(view, 'save')[0] as HTMLElement).click();
      await waitUntil(() => document.querySelector('.d4-dialog') !== null, 'Save dialog not shown', 5000);
      const saveBtn = Array.from(document.querySelectorAll('.d4-dialog button'))
        .find((b) => b.textContent?.trim() === 'Save') as HTMLElement;
      saveBtn.click();
      await waitUntil(() => savedRunId(view) != null, 'URL should carry the saved run id', 20000);
      try {
        expect((view.path ?? '').includes('currentStep='), false, 'currentStep must not be written in compact mode');
      } finally {
        await deleteSavedWorkflow(savedRunId(view)!);
      }
    });
  });
}

category('Editors: RFV saved runs', () => {
  const track = closeViewsAfter();

  test('Opening two runs quickly leaves the second one in the form', () => withRunStorage(async (f) => {
    await awaitWebComponents();
    const name = uniqueName('Quick open');
    // a run with an output table takes longer to load than one without
    const slow = await f.save(await f.dayRun({called: true}), {metadata: {name: `${name} slow`}});
    const fast = await f.saveDayRun(`${name} fast`);
    const view = track(await grok.functions.call('Compute2:RichFunctionViewEditor',
      {call: DG.Func.byName(DAY_CALCULATION).prepare()}) as unknown as DG.ViewBase);
    await awaitCheck(() => ribbonIcons(view, 'history').length > 0, 'Runs panel icon not shown', 15000);
    (ribbonIcons(view, 'history')[0] as HTMLElement).click();
    const row = (suffix: string) => Array.from(view.root.querySelectorAll<HTMLElement>('.c2-runs-list .c2-runs-row'))
      .find((el) => el.textContent?.includes(`${name} ${suffix}`));
    await awaitCheck(() => !!row('slow') && !!row('fast'), 'both runs should be listed', 15000);
    for (const suffix of ['slow', 'fast']) {
      for (const detail of [1, 2])
        row(suffix)!.dispatchEvent(new MouseEvent('click', {bubbles: true, detail}));
    }
    await awaitCheck(() => [slow.id, fast.id].includes(savedRunId(view) ?? ''), 'a run should open', 15000);
    // a load sent now comes back after the slower one that started before it
    await historyUtils.loadRun(slow.id);
    expect(savedRunId(view), fast.id);
  }));
});
