import * as grok from 'datagrok-api/grok';
// import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';


category('Scripts: Parameters', () => {
  const df = grok.data.demo.demog(20); 
  const col = df.getCol('sex');
  col.name = 'column';
  // const colList: string[] = ['age', 'weight'];
  const langs = ['Python', 'R', 'Octave', 'Grok', 'Julia', 'JavaScript', 'NodeJS'];

  for (const l of langs) {
    test(l, async () => {
      const f = await grok.functions.eval(`${_package.name}:${l}ParamsTest`);
      const call = f.prepare({i: 10, d: -20.1, b: false, s: 'abc', dt: '1961-06-01', df: df, col: col});
      await call.call();
      expect(call.getParamValue('ri'), 5);
      expect(call.getParamValue('rd'), 39.9);
      expect(call.getParamValue('rb'), true);
      expect(call.getParamValue('rs'), 'abcabc');
      expect(checkDate(call.getParamValue('rdt')), true);
      // expect(call.getParamValue('rdf').columns.length, 1);
    });
  }
});

function checkDate(dt: any): boolean {
  if (dt instanceof Date) return dt.toDateString() === 'Sun Jun 11 1961';
  return dt.format().includes('1961-06-11T00:00:00');
}
