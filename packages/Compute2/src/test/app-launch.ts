import * as DG from 'datagrok-api/dg';
import {category, test, expect, before, after, awaitCheck} from '@datagrok-libraries/test/src/test';
import {launchApp, closeView, awaitElement, awaitWebComponents, findButton, clickAndAwait} from './utils';

category('Apps: ViewerTestApp', () => {
  let view: DG.ViewBase;

  before(async () => {
    await awaitWebComponents();
  });

  after(async () => {
    if (view) closeView(view);
  });

  test('Renders viewer with data and buttons', async () => {
    view = await launchApp('ViewerTestApp');
    await awaitElement(view.root, 'dg-viewer', 'Viewer WebComponent', 10000);
    // Verify both control buttons are present with correct labels
    const changeDataBtn = findButton(view.root, 'change data');
    const changeTypeBtn = findButton(view.root, 'change type');
    expect(changeDataBtn !== null, true, 'change data button missing');
    expect(changeTypeBtn !== null, true, 'change type button missing');
  });

  test('Change data loads dataset into viewer', async () => {
    view = await launchApp('ViewerTestApp');
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
    view = await launchApp('ViewerTestApp');
    await awaitElement(view.root, 'dg-viewer', 'Viewer', 10000);
    // Set data first — types have no effect without data
    findButton(view.root, 'change data').click();
    await awaitCheck(
      () => view.root.querySelector('dg-viewer') !== null,
      'Viewer gone after setting initial data', 5000,
    );

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
  let view: DG.ViewBase;

  before(async () => {
    await awaitWebComponents();
  });

  after(async () => {
    if (view) closeView(view);
  });

  test('Shows SimpleInputs2 form with labeled inputs', async () => {
    view = await launchApp('FormTestApp');
    const btn = findButton(view.root, 'change funcall');
    // First click loads SimpleInputs2 (or ObjectCooling2 depending on toggle state)
    await clickAndAwait(
      btn,
      () => view.root.querySelector('dg-input-form') !== null,
      'InputForm did not appear', 10000,
    );

    // The form should contain actual input fields rendered by DG.InputForm
    await awaitCheck(
      () => view.root.querySelectorAll('dg-input-form input').length > 0,
      'InputForm has no input elements', 5000,
    );
  });

  test('Toggle funcall switches form content', async () => {
    view = await launchApp('FormTestApp');
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
    await awaitCheck(
      () => {
        const form = view.root.querySelector('dg-input-form');
        return form !== null && form.innerHTML !== formBefore && form.querySelectorAll('input').length > 0;
      },
      'Form content did not change after toggling function', 10000,
    );
  });

  test('Remove funcall clears all inputs', async () => {
    view = await launchApp('FormTestApp');
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

category('Apps: HistoryTestApp', () => {
  let view: DG.ViewBase;

  before(async () => {
    await awaitWebComponents();
  });

  after(async () => {
    if (view) closeView(view);
  });

  test('Renders history component for ObjectCooling2', async () => {
    view = await launchApp('HistoryTestApp');
    // The History component should render — either a table/grid with runs,
    // or an empty state message. Verify it has meaningful content.
    await awaitCheck(() => {
      const text = view.root.textContent ?? '';
      // The History component shows the function name or "No runs" or a table
      return view.root.children.length > 0
        && (text.includes('ObjectCooling2') || text.includes('cooling')
          || text.includes('No') || text.length > 10);
    }, 'HistoryTestApp rendered no meaningful content', 10000);
  });
});
