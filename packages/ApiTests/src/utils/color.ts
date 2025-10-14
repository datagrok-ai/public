// import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {
  intToHtmlA, intToHtml,
  intToRgba, intToRgb,
  htmlToInt, htmlToIntA,
} from '@datagrok-libraries/utils/src/color';
import * as DG from 'datagrok-api/dg';

category('Utils: color', () => {
  const data: { [testName: string]: { src: any, tgt: any, to: (src: any) => any } } = {
    intToHex: {src: 0x00bfff, tgt: '#00bfff', to: intToHtml},
    intToHexA: {src: 0x3300bfff, tgt: '#00bfff33', to: intToHtmlA},
    intToRgb: {src: 0x00bfff, tgt: 'rgb(0, 191, 255)', to: intToRgb},
    intToRgbA: {src: 0x3300bfff, tgt: `rgba(0, 191, 255, ${0x33 / 255})`, to: intToRgba},
    htmlToInt: {src: '#00bfff', tgt: 0x00bfff, to: htmlToInt},
    htmlToIntA: {src: '#00bfff33', tgt: 0x3300bfff, to: htmlToIntA},
  };

  for (const [testName, testCfg] of Object.entries(data)) {
    test(testName, async () => {
      const res = testCfg.to(testCfg.src);
      expect(testCfg.tgt, res);
    });
  }
  
  test('DG.toHtml', async () => {
    const res = DG.Color.toHtml(10);
    expect("#00000a", res);
  });
});
