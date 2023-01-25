import {before, category, test, expect} from '@datagrok-libraries/utils/src/test';
import { ColumnList } from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { _package } from '../package-test';
import * as DG from 'datagrok-api/dg';

category('Scripts Params', () => {

  const df = grok.data.demo.demog(20); 
  const col = df.col('sex');
  let colList:string[] = ['age', 'weight'];

  test('python.test', async () => {
    let f = await grok.functions.eval(`${_package.name}:PythonParamsTest`);

    let call = f.prepare({i:10, d: -20.1, b: false, s: 'abc', dt:'1992-09-20 00:00:00', df:df, col:col});
    await call.call();
    
    expect(call.getParamValue('ri'), 5);
    expect(call.getParamValue('rd'), 39.9);
    expect(call.getParamValue('rb'), true);
    expect(call.getParamValue('rs'), 'abcabc');
    expect(call.getParamValue('rdt').toString(),'Wed, 09 Sep 1992 21:00:00 GMT');
    expect(call.getParamValue('rdf').columns.length, 1);    
  });

  test('r.test', async () => {
    let f = await grok.functions.eval(`${_package.name}:RParamsTest`);

    let call = f.prepare({i:10, d: -20.1, b: false, s: 'abc', dt:'1992-09-20 00:00:00', df:df, col:col});
    await call.call();
    
    expect(call.getParamValue('ri'), 5);
    expect(call.getParamValue('rd'), 39.9);
    expect(call.getParamValue('rb'), true);
    expect(call.getParamValue('rs'), 'abcabc');
    expect(call.getParamValue('rdt').toString(),'Wed, 09 Sep 1992 21:00:00 GMT');
    expect(call.getParamValue('rdf').columns.length, 1);    
  });

  test('octave.test', async () => {
    let f = await grok.functions.eval(`${_package.name}:OctaveParamsTest`);

    let call = f.prepare({i:10, d: -20.1, b: false, s: 'abc', df:df, col:col, col_list: colList});
    await call.call();
    
    expect(call.getParamValue('ri'), 73);
    expect(call.getParamValue('rd'), 9.9);
    expect(call.getParamValue('rb'), true);
    expect(call.getParamValue('rs'), 'abc-age');
    expect(call.getParamValue('rdf').columns.length, 2); 
  });
});
