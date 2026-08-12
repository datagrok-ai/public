import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

// Pins the Dart behavior the dg-viewer webcomponent relies on: a frame swap followed
// synchronously by consistent options never evaluates a torn look/frame pair, because
// the frame swap's look refresh is deferred while setOptions corrects the look first.
// If either direction starts failing here, the platform made the frame-swap refresh
// eager and the webcomponent's ViewerHost swap needs a redesign.

const makeDf = (withSpecies: boolean) => {
  const columns = [
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'time', [0, 1, 2, 3]),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'value', [10, 20, 30, 40]),
    ...withSpecies ?
      [DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'species', ['a', 'a', 'b', 'b'])] : [],
  ];
  return DG.DataFrame.fromColumns(columns);
};

const captureNullErrors = async (action: () => void): Promise<string[]> => {
  const errors: string[] = [];
  const original = console.error;
  console.error = (...args: any[]) => {
    const message = args.map((a) => `${a}`).join(' ');
    if (message.includes('NullError') || message.includes('method not found') ||
      message.includes('Cannot read properties of null'))
      errors.push(message);
    original.apply(console, args);
  };
  try {
    action();
    await new Promise((resolve) => setTimeout(resolve, 400));
  } finally {
    console.error = original;
  }
  return errors;
};

category('RunComparison: viewer look swap pins', () => {
  test('split deselect: frame swap then split-free options is clean', async () => {
    const view = grok.shell.newView('swap-pin-deselect');
    try {
      const viewer = await makeDf(true).plot.fromType(DG.VIEWER.LINE_CHART) as DG.Viewer;
      viewer.setOptions({xColumnName: 'time', yColumnNames: ['value'], splitColumnNames: ['species']});
      view.append(viewer.root);
      await new Promise((resolve) => setTimeout(resolve, 300));
      const errors = await captureNullErrors(() => {
        viewer.dataFrame = makeDf(false);
        viewer.setOptions({xColumnName: 'time', yColumnNames: ['value']});
      });
      expect(errors.length, 0);
    } finally {
      view.close();
    }
  });

  test('split select: frame swap then split options is clean', async () => {
    const view = grok.shell.newView('swap-pin-select');
    try {
      const viewer = await makeDf(false).plot.fromType(DG.VIEWER.LINE_CHART) as DG.Viewer;
      viewer.setOptions({xColumnName: 'time', yColumnNames: ['value']});
      view.append(viewer.root);
      await new Promise((resolve) => setTimeout(resolve, 300));
      const errors = await captureNullErrors(() => {
        viewer.dataFrame = makeDf(true);
        viewer.setOptions({xColumnName: 'time', yColumnNames: ['value'], splitColumnNames: ['species']});
      });
      expect(errors.length, 0);
    } finally {
      view.close();
    }
  });
});
