import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getWrapSettings, setMonomersPerRow} from '@datagrok-libraries/hwe';

import {initHelmMainPackage} from './utils';
import {getWrapWidthWidget} from '../widgets/wrap-width-widget';
import {resolveMonomersPerRow} from '../utils/wrap-width';
import {DEFAULT_MONOMERS_PER_ROW} from '../constants';

// The int input inside the widget produced by `getWrapWidthWidget`.
function widgetInput(widget: DG.Widget): HTMLInputElement {
  const input = widget.root.querySelector('input') as HTMLInputElement | null;
  if (input === null) throw new Error('wrap-width widget: no input element');
  return input;
}

// Set the input the way a user would, so the widget's change handler runs.
function typeValue(input: HTMLInputElement, value: number): void {
  input.value = String(value);
  input.dispatchEvent(new Event('input', {bubbles: true}));
  input.dispatchEvent(new Event('change', {bubbles: true}));
}

category('WrapWidth: package property', () => {
  test('reads an int property', async () => {
    expect(resolveMonomersPerRow(30), 30);
  });

  test('reads a string defaultValue straight from package.json', async () => {
    // An `int` property that has not round-tripped through the server arrives
    // as the string written in package.json.
    expect(resolveMonomersPerRow('30'), 30);
    expect(resolveMonomersPerRow('20'), 20);
  });

  test('falls back when settings are missing or unusable', async () => {
    // A NaN reaching the layout would produce a drawing with no coordinates,
    // so every unusable value must land on the default instead.
    for (const raw of [undefined, null, '', 'abc', Number.NaN, 0, -5, {}])
      expect(resolveMonomersPerRow(raw), DEFAULT_MONOMERS_PER_ROW);
  });

  test('default matches the value declared in package.json', async () => {
    expect(DEFAULT_MONOMERS_PER_ROW, 20);
  });
});

category('WrapWidth: session widget', () => {
  // The wrap width is process-wide state inside hwe; put it back afterwards so
  // a failure here cannot change how the rest of the suite renders.
  let restore: number;

  before(async () => {
    await initHelmMainPackage();
    restore = getWrapSettings().monomersPerRow;
  });

  after(async () => {
    setMonomersPerRow(restore);
  });

  test('opens showing the current session value', async () => {
    setMonomersPerRow(14);
    const widget = getWrapWidthWidget();
    try {
      expect(Number(widgetInput(widget).value), 14);
    } finally {
      widget.detach();
    }
  });

  test('editing the input changes the session value', async () => {
    setMonomersPerRow(DEFAULT_MONOMERS_PER_ROW);
    const widget = getWrapWidthWidget();
    try {
      typeValue(widgetInput(widget), 8);
      expect(getWrapSettings().monomersPerRow, 8);
    } finally {
      widget.detach();
    }
  });

  test('editing does NOT change the package default', async () => {
    // This is the whole point of the split: the package property stays the
    // value every new session starts from.
    const widget = getWrapWidthWidget();
    try {
      typeValue(widgetInput(widget), 5);
      expect(getWrapSettings().monomersPerRow, 5);
      expect(resolveMonomersPerRow('20'), 20);
      expect(DEFAULT_MONOMERS_PER_ROW, 20);
    } finally {
      widget.detach();
    }
  });

  test('reflects a change made elsewhere', async () => {
    const widget = getWrapWidthWidget();
    try {
      setMonomersPerRow(13);
      expect(Number(widgetInput(widget).value), 13);
    } finally {
      widget.detach();
    }
  });

  test('stops tracking after detach', async () => {
    const widget = getWrapWidthWidget();
    const input = widgetInput(widget);
    typeValue(input, 9);
    widget.detach();
    setMonomersPerRow(17);
    // The subscription goes through `widget.sub`, so `detach` unsubscribes it.
    expect(Number(input.value), 9);
    expect(getWrapSettings().monomersPerRow, 17);
  });

  test('a rejected value leaves the session value alone', async () => {
    setMonomersPerRow(11);
    const widget = getWrapWidthWidget();
    try {
      const input = widgetInput(widget);
      input.value = '';
      input.dispatchEvent(new Event('input', {bubbles: true}));
      input.dispatchEvent(new Event('change', {bubbles: true}));
      expect(getWrapSettings().monomersPerRow, 11);
    } finally {
      widget.detach();
    }
  });
});
