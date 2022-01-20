import { after, before, category, delay, expect, test } from "@datagrok-libraries/utils/src/test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { checkHTMLElement } from './utils';

category('UI: Range slider', () => {
  let v: DG.View;

  let minRange = 0;
  let maxRange = 10;
  let min = 2;
  let max = 5;
  let range = ui.rangeSlider(minRange, maxRange, min, max);

  before(async () => {
    v = grok.shell.newView('');
  });

  test('range.minRange', async () => {
    expect(range.minRange, minRange)
  })

  test('range.maxRange', async () => {
    expect(range.maxRange, maxRange)
  })

  test('range.min', async () => {
    expect(range.min, min)
  })

  test('range.max', async () => {
    expect(range.max, max)
  })

  test('range.root', async () => {
    checkHTMLElement('range', range.root, v, '.d4-range-selector');
  });

  test('range.onValuesChanged', async () => {
    let check = false;
    try {
      range.onValuesChanged.subscribe((_) => { check = true });
      range.setValues(1, 9, 3, 4);
      if (check = false)
        throw 'onValuesChanged error';
    }
    catch (x) {
      throw 'rangeSlider: ' + x
    }
    finally {
      range.setValues(minRange, maxRange, min, max);
    }
  })

  after(async () => {
    v.close();
    grok.shell.closeAll();
  });

});