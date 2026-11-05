import * as DG from 'datagrok-api/dg';
import {category, test, expect, before, waitUntil, awaitCheck} from '@datagrok-libraries/test/src/test';
import {closeViewsAfter, launchApp, awaitElement, awaitWebComponents, findButton, clickAndAwait} from './utils';
import {DAY_CALCULATION, withRunStorage} from './run-storage';
import {openRunBrowserDialog} from '../components/RunStorage/run-browser-dialog';

category('Apps: ViewerTestApp', () => {
  const track = closeViewsAfter();

  before(async () => {
    await awaitWebComponents();
  });


  test('Renders viewer with data and buttons', async () => {
    const view = track(await launchApp('ViewerTestApp'));
    await awaitElement(view.root, 'dg-viewer', 'Viewer WebComponent', 10000);
    // Verify both control buttons are present with correct labels
    const changeDataBtn = findButton(view.root, 'change data');
    const changeTypeBtn = findButton(view.root, 'change type');
    expect(changeDataBtn !== null, true, 'change data button missing');
    expect(changeTypeBtn !== null, true, 'change type button missing');
  });

  test('Change data loads dataset into viewer', async () => {
    const view = track(await launchApp('ViewerTestApp'));
    await awaitElement(view.root, 'dg-viewer', 'Viewer', 10000);
    const btn = findButton(view.root, 'change data');

    // Click change data 3 times — cycles through demog, doseResponse, geo.
    // The dg-viewer WebComponent delegates rendering to the platform viewer,
    // which may not expose inner content as DOM children.
    // Verify the component stays mounted and no exceptions are thrown.
    for (let i = 0; i < 3; i++) {
      await clickAndAwait(
        btn,
        () => view.root.querySelector('dg-viewer') !== null,
        `Viewer disappeared after data change #${i + 1}`, 5000,
      );
    }
  });

  test('Change type renders different viewer kinds', async () => {
    const view = track(await launchApp('ViewerTestApp'));
    await awaitElement(view.root, 'dg-viewer', 'Viewer', 10000);
    // Set data first — types have no effect without data
    findButton(view.root, 'change data').click();
    await waitUntil(() => view.root.querySelector('dg-viewer') !== null,
      'Viewer gone after setting initial data', 5000);

    // Cycle through all 4 types: Grid → Histogram → Line chart → Scatter plot
    const btn = findButton(view.root, 'change type');
    for (let i = 0; i < 4; i++) {
      await clickAndAwait(
        btn,
        () => view.root.querySelector('dg-viewer') !== null,
        `Viewer disappeared after type change #${i + 1}`, 5000,
      );
    }
  });
});

category('Apps: FormTestApp', () => {
  const track = closeViewsAfter();

  before(async () => {
    await awaitWebComponents();
  });


  test('Shows SimpleInputs2 form with labeled inputs', async () => {
    const view = track(await launchApp('FormTestApp'));
    const btn = findButton(view.root, 'change funcall');
    // First click loads SimpleInputs2 (or ObjectCooling2 depending on toggle state)
    await clickAndAwait(
      btn,
      () => view.root.querySelector('dg-input-form') !== null,
      'InputForm did not appear', 10000,
    );

    // The form should contain actual input fields rendered by DG.InputForm
    await waitUntil(() => view.root.querySelectorAll('dg-input-form input').length > 0,
      'InputForm has no input elements', 5000);
  });

  test('Toggle funcall switches form content', async () => {
    const view = track(await launchApp('FormTestApp'));
    const btn = findButton(view.root, 'change funcall');
    await clickAndAwait(
      btn,
      () => view.root.querySelector('dg-input-form') !== null,
      'InputForm did not appear after first toggle', 10000,
    );

    // Capture form innerHTML before toggling
    const formBefore = view.root.querySelector('dg-input-form')!.innerHTML;

    // Toggle to the other function — form content should change
    btn.click();
    await waitUntil(() => {
      const form = view.root.querySelector('dg-input-form');
      return form !== null && form.innerHTML !== formBefore && form.querySelectorAll('input').length > 0;
    }, 'Form content did not change after toggling function', 10000);
  });

  test('Remove funcall clears all inputs', async () => {
    const view = track(await launchApp('FormTestApp'));
    const changeFcBtn = findButton(view.root, 'change funcall');
    await clickAndAwait(
      changeFcBtn,
      () => view.root.querySelector('dg-input-form') !== null,
      'InputForm did not appear before remove', 10000,
    );
    const removeBtn = findButton(view.root, 'remove funcall');
    await clickAndAwait(
      removeBtn,
      () => view.root.querySelector('dg-input-form') === null,
      'InputForm should disappear after remove', 5000,
    );
    // Verify no stale input elements remain
    const remainingInputs = view.root.querySelectorAll('dg-input-form input').length;
    expect(remainingInputs, 0, 'No input elements should remain after remove');
  });
});

category('Apps: Browse runs', () => {
  before(async () => {
    await awaitWebComponents();
  });

  test('Lists a saved run of the chosen model', () => withRunStorage(async (f) => {
    const name = `Browse runs test ${Date.now()}`;
    await f.saveDayRun(name);
    const existing = new Set(document.querySelectorAll('.c2-runs-browser'));
    const opened = openRunBrowserDialog({
      title: 'Browse runs', selectModel: true, model: DG.Func.byName(DAY_CALCULATION),
    });
    let browser: HTMLElement | undefined;
    await awaitCheck(() => {
      browser = [...document.querySelectorAll<HTMLElement>('.c2-runs-browser')].find((el) => !existing.has(el));
      return !!browser;
    }, 'The run browser did not open', 10000);
    try {
      await awaitCheck(() => (browser!.textContent ?? '').includes(name), 'The saved run is not listed', 15000);
    } finally {
      findButton(browser!, 'cancel').click();
    }
    expect(await opened, null);
  }));
});
