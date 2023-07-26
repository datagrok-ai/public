import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';

category(`Scripts`, () => {
  test('Outputs conversion', async () => {
    const out = await grok.functions.call('CVMTests:DummyPython');
    expect(out.df1 instanceof DG.DataFrame, true); // works with DG.toJs(out.df1)
    expect(out.df2 instanceof DG.DataFrame, true); // works with DG.toJs(out.df2)
  });
});

const langs = ['Python', 'R', 'Octave', 'Grok', 'Julia', 'JavaScript']; // NodeJS, Julia

for (const l of langs) {
  category(`Script parameters: ${l}`, () => {
    let call: DG.FuncCall;

    before(async () => {
      const df = grok.data.demo.demog(20);
      const f = await grok.functions.eval(`${_package.name}:${l}ParamsTest`);
      const params = {i: 10, d: -20.1, b: false, s: 'abc', df: df};
      if (l !== 'Octave') Object.assign(params, {dt: '1961-06-01'}); // GROK-12417
      if (!['R', 'Julia', 'NodeJS', 'Octave'].includes(l)) Object.assign(params, {m: {a: 1}}); // GROK-12452
      call = f.prepare(params);
      await call.call();
    });

    test('int', async () => {
      if (call.getParamValue('ri') !== 5) throw new Error(`${call.getParamValue('ri')} != 5`);
    });
    test('double', async () => {
      if (call.getParamValue('rd') !== 39.9) throw new Error(`${call.getParamValue('rd')} != 39.9`);
    }, {skipReason: l === 'Julia' ? 'GROK-12386' : undefined});
    test('bool', async () => {
      if (!call.getParamValue('rb')) throw new Error(`${call.getParamValue('rb')} != true`);
    }, {skipReason: l === 'Julia' ? 'GROK-12386' : undefined});
    test('string', async () => {
      if (call.getParamValue('rs') !== 'abcabc') throw new Error(`"${call.getParamValue('rs')}" != "abcabc"`);
    }, {skipReason: l === 'Julia' ? 'GROK-12386' : undefined});
    test('datetime', async () => {
      if (!checkDate(call.getParamValue('rdt')))
        throw new Error(`${call.getParamValue('rdt')} != 1961-06-11`);
    }, {skipReason: ['Octave', 'Python', 'Julia'].includes(l) ? 'GROK-12417' : undefined});
    test('map', async () => {
      if (call.getParamValue('rm').b !== 5) throw new Error(`${call.getParamValue('rm').b} != 5`);
    }, {skipReason: ['R', 'Julia', 'NodeJS', 'Octave'].includes(l) ? 'GROK-12452' : undefined});
    test('dataframe', async () => {
      if (call.getParamValue('rdf').columns.length !== 11 || call.getParamValue('rdf').rowCount !== 20) {
        throw new Error(`df shape (${call.getParamValue('rdf').columns.length}, ${
          call.getParamValue('rdf').rowCount}) != (11, 20)`);
      }
    }, {skipReason: l === 'Julia' ? 'GROK-12386' : undefined});
  });
}

function checkDate(dt: any): boolean {
  if (dt == null) return false;
  if (dt instanceof Date) return dt.toDateString() === 'Sun Jun 11 1961';
  return dt.format().includes('1961-06-11T');
}
