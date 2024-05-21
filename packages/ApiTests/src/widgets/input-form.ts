import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {before, category, expect, test, after, assure, expectArray, expectTable} from '@datagrok-libraries/utils/src/test';

category('ValueLookup with no nullables', () => {
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let form: DG.InputForm;

  before(async () => {
    funcCall = (await grok.functions.eval('ApiTests:ValueLookup')).prepare();
  
    form = await DG.InputForm.forFuncCall(funcCall);
    inputs = {
      'model': form.getInput('model') ?? null,
      'mpg': form.getInput('mpg') ?? null,
      'cyl': form.getInput('cyl') ?? null,
      'disp': form.getInput('disp') ?? null,
    };
  });
  
  test('inputform created', async () => {
    assure.notNull(form);
  });

  test('lookup items', async () => {
    expectArray((inputs['model'] as DG.ChoiceInput<string>).items, 
      ['Mazda RX4', 'Mazda RX4 Wag', 'Datsun 710', 'Hornet 4 Drive', 'Hornet Sportabout']);
  });
  
  test('initial values', async () => {
    expect(inputs['model'].value, 'Mazda RX4');
    expect(inputs['mpg'].value, 21);
    expect(inputs['cyl'].value, 6);
    expect(inputs['disp'].value, 160);
  });
});

category('ValueLookup with with nullables', () => {
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let form: DG.InputForm;

  before(async () => {
    funcCall = (await grok.functions.eval('ApiTests:ValueLookupWithNullable')).prepare();
  
    form = await DG.InputForm.forFuncCall(funcCall);
    inputs = {
      'model': form.getInput('model') ?? null,
      'mpg': form.getInput('mpg') ?? null,
      'cyl': form.getInput('cyl') ?? null,
      'disp': form.getInput('disp') ?? null,
    };
  });
  
  test('inputform created', async () => {
    assure.notNull(form);
  });

  test('lookup items', async () => {
    expectArray((inputs['model'] as DG.ChoiceInput<string>).items, 
      ['', 'Mazda RX4', 'Mazda RX4 Wag', 'Datsun 710', 'Hornet 4 Drive', 'Hornet Sportabout']);
  });
  
  test('initial values', async () => {
    expect(inputs['model'].value, null);
    expect(inputs['mpg'].value, null);
    expect(inputs['cyl'].value, null);
    expect(inputs['disp'].value, null);
  });
});

category('InputForm API', () => {
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let newFuncCall: DG.FuncCall;
  let form: DG.InputForm;

  const demog = grok.data.demo.demog(10);

  before(async () => {
    funcCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
  
    form = await DG.InputForm.forFuncCall(funcCall, {twoWayBinding: true});
    inputs = {
      'stringInput': form.getInput('stringInput') ?? null,
      'intInput': form.getInput('intInput') ?? null,
      'doubleInput': form.getInput('doubleInput') ?? null,
      'boolInput': form.getInput('boolInput') ?? null,
      'choiceInput': form.getInput('choiceInput') ?? null,
      'tableInput': form.getInput('tableInput') ?? null,
    };
  });
  
  test('form to funccall bind', async () => {
    inputs['stringInput'].value = 'test2';
    inputs['intInput'].value = 4;
    inputs['doubleInput'].value = 4.14;
    inputs['boolInput'].value = false;
    inputs['choiceInput'].value = '2';
    inputs['tableInput'].value = demog;

    expect(funcCall.inputs['stringInput'], 'test2');
    expect(funcCall.inputs['intInput'], 4);
    expect(funcCall.inputs['doubleInput'], 4.14);
    expect(funcCall.inputs['boolInput'], false);
    expect(funcCall.inputs['choiceInput'], '2');
    expectTable(funcCall.inputs['tableInput'], demog);
  });
  
  test('funcall to form bind', async () => {
    const geo = grok.data.demo.geo(10);

    funcCall.inputs['stringInput'] = 'test';
    funcCall.inputs['intInput'] = 3;
    funcCall.inputs['doubleInput'] = 3.14;
    funcCall.inputs['boolInput'] = true;
    funcCall.inputs['choiceInput'] = '1';
    funcCall.inputs['tableInput'] = geo;

    expect(inputs['stringInput'].value, 'test');
    expect(inputs['intInput'].value, 3);
    expect(inputs['doubleInput'].value, 3.14);
    expect(inputs['boolInput'].value, true);
    expect(inputs['choiceInput'].value, '1');
    expectTable(inputs['tableInput'].value, geo);
  });

  test('form on input change observable', async () => {
    const changedInputPropNames = [] as string[];
    const changeSub = form.onInputChanged.subscribe((propName) => {
      changedInputPropNames.push(propName);
    });

    inputs['stringInput'].value = 'test2';
    inputs['intInput'].value = 4;
    inputs['doubleInput'].value = 4.14;
    inputs['boolInput'].value = false;
    inputs['choiceInput'].value = '2';
    inputs['tableInput'].value = demog;

    expectArray(
      changedInputPropNames, 
      ['stringInput', 'intInput', 'doubleInput', 'boolInput', 'choiceInput', 'tableInput'],
    );
    changeSub.unsubscribe();
  });

  test('source funccall replacement', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;

    expect(inputs['stringInput'].value, 'test');
    expect(inputs['intInput'].value, 3);
    expect(inputs['doubleInput'].value, 3.14);
    expect(inputs['boolInput'].value, true);
    expect(inputs['choiceInput'].value, '1');
    expect(inputs['tableInput'].value, null); // no default value here
  });

  test('form to funccall bind after replace', async () => {
    inputs['stringInput'].value = 'test2';
    inputs['intInput'].value = 4;
    inputs['doubleInput'].value = 4.14;
    inputs['boolInput'].value = false;
    inputs['choiceInput'].value = '2';
    inputs['tableInput'].value = demog;

    expect(newFuncCall.inputs['stringInput'], 'test2');
    expect(newFuncCall.inputs['intInput'], 4);
    expect(newFuncCall.inputs['doubleInput'], 4.14);
    expect(newFuncCall.inputs['boolInput'], false);
    expect(newFuncCall.inputs['choiceInput'], '2');
    expectTable(newFuncCall.inputs['tableInput'], demog);
  });
  
  test('funccall to form bind after replace', async () => {
    const geo = grok.data.demo.geo(10);

    newFuncCall.inputs['stringInput'] = 'test';
    newFuncCall.inputs['intInput'] = 3;
    newFuncCall.inputs['doubleInput'] = 3.14;
    newFuncCall.inputs['boolInput'] = true;
    newFuncCall.inputs['choiceInput'] = '1';
    newFuncCall.inputs['tableInput'] = geo;

    expect(inputs['stringInput'].value, 'test');
    expect(inputs['intInput'].value, 3);
    expect(inputs['doubleInput'].value, 3.14);
    expect(inputs['boolInput'].value, true);
    expect(inputs['choiceInput'].value, '1');
    expectTable(inputs['tableInput'].value, geo);
  });

  test('form on input change observable after replace', async () => {
    const changedInputPropNames = [] as string[];
    const changeSub = form.onInputChanged.subscribe((propName) => {
      changedInputPropNames.push(propName);
    });

    inputs['stringInput'].value = 'test2';
    inputs['intInput'].value = 4;
    inputs['doubleInput'].value = 4.14;
    inputs['boolInput'].value = false;
    inputs['choiceInput'].value = '2';
    inputs['tableInput'].value = demog;

    expectArray(
      changedInputPropNames, 
      ['stringInput', 'intInput', 'doubleInput', 'boolInput', 'choiceInput', 'tableInput'],
    );
    changeSub.unsubscribe();
  });
});
