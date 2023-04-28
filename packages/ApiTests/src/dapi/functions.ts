import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

const GDF = grok.dapi.functions;

category('Dapi: functions.calls', async () => {
  const xValue = 1.5;

  test('save', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    funcCall.newId();
    const savedFuncCall = await GDF.calls.save(funcCall);
    expect(savedFuncCall.inputs['x'], funcCall.inputs['x']);
  });

  test('save options', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    funcCall.newId();
    funcCall.options['testName'] = 'testValue';
    const savedFuncCall = await GDF.calls.save(funcCall);
    expect(savedFuncCall.options['testName'], 'testValue');
  });

  test('save parentFunccall', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    const parentFunc: DG.Func = await grok.functions.eval('Cos');
    const parentFuncCall = await parentFunc.prepare({x: xValue}).call();
    parentFuncCall.newId();
    funcCall.newId();
    funcCall.parentCall = parentFuncCall;
    await GDF.calls.save(parentFuncCall);
    await GDF.calls.save(funcCall);

    const loadedFunCall = await GDF.calls.include('parentCall').find(funcCall.id);
    expect(loadedFunCall.parentCall.id, parentFuncCall.id);
  });

  test('load package function call', async () => {
    const packFunc: DG.Func = await grok.functions.eval('ApiTests:dummyPackageFunction');
    const funcCall = await packFunc.prepare({a: 1, b: 2}).call();
    funcCall.newId();
    await GDF.calls.save(funcCall);

    // expect no-throw
    await GDF.calls.find(funcCall.id);
  });

  test('load script call', async () => {
    const packFunc: DG.Func = await grok.functions.eval('ApiTests:dummyPackageScript');
    const funcCall = await packFunc.prepare({a: 1, b: 2}).call();
    funcCall.newId();
    await GDF.calls.save(funcCall);

    // expect no-throw
    await GDF.calls.find(funcCall.id);
  });

  test('list', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    funcCall.newId();
    await GDF.calls.save(funcCall);
    const loadedFuncCalls = await GDF.calls.filter(`func.id="${funcCall.func.id}"`).list();
    expect(loadedFuncCalls.some((loadedCall) => loadedCall.id === funcCall.id), true);
  });

  test('find', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    funcCall.newId();
    await GDF.calls.save(funcCall);
    const loadedFuncCall = await GDF.calls.find(funcCall.id);
    expect(loadedFuncCall.inputs['x'], xValue);
  });

  test('find options', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    funcCall.newId();
    funcCall.options['testName'] = 'testValue';
    await GDF.calls.save(funcCall);
    const loadedFuncCall = await GDF.calls.include('options').find(funcCall.id);
    expect(loadedFuncCall.options['testName'], 'testValue');
  });

  test('delete', async () => {
    const func: DG.Func = await grok.functions.eval('Sin');
    const funcCall = await func.prepare({x: xValue}).call();
    expect((await GDF.calls.filter(`func.id="${funcCall.func.id}"`).list()).includes(funcCall), false);
    await GDF.calls.delete(funcCall);
    expect((await GDF.calls.filter(`func.id="${funcCall.func.id}"`).list()).includes(funcCall), false);
  });
});
