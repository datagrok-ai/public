import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test, assure, expectArray, expectTable} from '@datagrok-libraries/utils/src/test';
import {take} from 'rxjs/operators';

const zakievAufar = 'aufar.zakiev@softwarecountry.com';
const createForm = async (withNullable: boolean = false) => {
  const funcCall: DG.FuncCall = (await grok.functions.eval(`ApiTests:ValueLookup${withNullable ? 'WithNullable' : ''}`))
    .prepare();
  const form: DG.InputForm = await DG.InputForm.forFuncCall(funcCall);
  const inputs: Record<string, DG.InputBase> = {
    'model': form.getInput('model') ?? null,
    'mpg': form.getInput('mpg') ?? null,
    'cyl': form.getInput('cyl') ?? null,
    'disp': form.getInput('disp') ?? null,
  };
  return {form, inputs};
};

category('Widgets: ValueLookup with no nullables', () => {
  test('inputform created', async () => {
    assure.notNull(await createForm());
  });

  test('lookup items', async () => {
    const {inputs} = await createForm();
    expectArray((inputs['model'] as DG.ChoiceInput<string>).items,
      ['Mazda RX4', 'Mazda RX4 Wag', 'Datsun 710', 'Hornet 4 Drive', 'Hornet Sportabout']);
  });

  test('initial values', async () => {
    const {inputs} = await createForm();
    expect(inputs['model'].value, 'Mazda RX4');
    expect(inputs['mpg'].value, 21);
    expect(inputs['cyl'].value, 6);
    expect(inputs['disp'].value, 160);
  });
}, {owner: 'dkovalyov@datagrok.ai'});

category('Widgets: ValueLookup with with nullables', () => {
  test('inputform created', async () => {
    assure.notNull(await createForm(true));
  });

  test('lookup items', async () => {
    const {inputs} = await createForm(true);
    expectArray((inputs['model'] as DG.ChoiceInput<string>).items,
      [null, 'Mazda RX4', 'Mazda RX4 Wag', 'Datsun 710', 'Hornet 4 Drive', 'Hornet Sportabout']);
  });

  test('initial values', async () => {
    const {inputs} = await createForm(true);
    expect(inputs['model'].value, null);
    expect(inputs['mpg'].value, null);
    expect(inputs['cyl'].value, null);
    expect(inputs['disp'].value, null);
  });
}, {owner: 'dkovalyov@datagrok.ai'});

category('Widgets: InputForm fc replacement edge cases', () => {
  test('inputform created', async () => {
    assure.notNull(await createForm());
  });

  test('initial values', async () => {
    const {inputs} = await createForm();
    expect(inputs['model'].value, 'Mazda RX4');
    expect(inputs['mpg'].value, 21);
    expect(inputs['cyl'].value, 6);
    expect(inputs['disp'].value, 160);
  }, {owner: zakievAufar});

  test('source replace w/o value lookup run', async () => {
    let {form, inputs} = await createForm();
    const newFuncCall = (await grok.functions.eval('ApiTests:ValueLookup')).prepare({
      model: 'Mazda RX4',
      with_choices: '0',
    });
    form.source = newFuncCall;
    inputs = {
      'model': form.getInput('model') ?? null,
      'mpg': form.getInput('mpg') ?? null,
      'cyl': form.getInput('cyl') ?? null,
      'disp': form.getInput('disp') ?? null,
    };
    expect(inputs['model'].value, 'Mazda RX4');
    expect(inputs['mpg'].value, null);
    expect(inputs['cyl'].value, null);
    expect(inputs['disp'].value, null);
    expect(inputs['with_choices'].value, '0');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15741', owner: zakievAufar});
}, {owner: 'dkovalyov@datagrok.ai'});

category('Widgets: InputForm API', () => {
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let newFuncCall: DG.FuncCall;
  let form: DG.InputForm;

  const demog = grok.data.demo.demog(10);

  function updateInputs(): void {
    inputs = {
      'stringInput': form.getInput('stringInput') ?? null,
      'intInput': form.getInput('intInput') ?? null,
      'doubleInput': form.getInput('doubleInput') ?? null,
      'boolInput': form.getInput('boolInput') ?? null,
      'choiceInput': form.getInput('choiceInput') ?? null,
      'tableInput': form.getInput('tableInput') ?? null,
    };
  }

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
  }, {owner: zakievAufar});

  test('funcall to form bind', async () => {
    const geo = grok.data.demo.geo(10);

    funcCall.inputs['stringInput'] = 'test2';
    funcCall.inputs['intInput'] = 4;
    funcCall.inputs['doubleInput'] = 4.14;
    funcCall.inputs['boolInput'] = false;
    funcCall.inputs['choiceInput'] = '2';
    funcCall.inputs['tableInput'] = geo;

    expect(inputs['stringInput'].value, 'test2');
    expect(inputs['intInput'].value, 4);
    expect(inputs['doubleInput'].value, 4.14);
    expect(inputs['boolInput'].value, false);
    expect(inputs['choiceInput'].value, '2');
    expectTable(inputs['tableInput'].value, geo);
  }, {owner: zakievAufar});

  test('form on input change observable', async () => {
    const changedInputPropNames = [] as string[];
    form.onInputChanged.pipe(take(1)).subscribe((ed) => {
      changedInputPropNames.push(ed.args.input.property.name);
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
  }, {owner: zakievAufar});

  test('source funccall replacement', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare({'stringInput': 'test2'});
    form.source = newFuncCall;
    updateInputs();

    expect(inputs['stringInput'].value, 'test2');
  }, {owner: zakievAufar});

  test('form to funccall bind after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;
    updateInputs();

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
  }, {owner: zakievAufar});

  test('funccall to form bind after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;

    const geo = grok.data.demo.geo(10);

    newFuncCall.inputs['stringInput'] = 'test';
    newFuncCall.inputs['intInput'] = 3;
    newFuncCall.inputs['doubleInput'] = 3.14;
    newFuncCall.inputs['boolInput'] = true;
    newFuncCall.inputs['choiceInput'] = '1';
    newFuncCall.inputs['tableInput'] = geo;
    updateInputs();

    expect(inputs['stringInput'].value, 'test');
    expect(inputs['intInput'].value, 3);
    expect(inputs['doubleInput'].value, 3.14);
    expect(inputs['boolInput'].value, true);
    expect(inputs['choiceInput'].value, '1');
    expectTable(inputs['tableInput'].value, geo);
  }, {skipReason: 'GROK-16408', owner: zakievAufar});

  test('form on input change observable after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;
    updateInputs();

    const changedInputPropNames = [] as string[];
    form.onInputChanged.pipe(take(1)).subscribe((ed) => {
      changedInputPropNames.push(ed.args.input.property.name);
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
  }, {owner: zakievAufar});
}, {owner: 'dkovalyov@datagrok.ai'});


category('Widgets: InputForm w/ custom input', () => {
  let inputs = {} as Record<string, DG.InputBase>;
  let funcCall: DG.FuncCall;
  let newFuncCall: DG.FuncCall;
  let form: DG.InputForm;

  
  function updateInputs(): void {
    inputs = {
      'stringInput': form.getInput('stringInput') ?? null,
      'intInput': form.getInput('intInput') ?? null,
      'doubleInput': form.getInput('doubleInput') ?? null,
      'boolInput': form.getInput('boolInput') ?? null,
      'choiceInput': form.getInput('choiceInput') ?? null,
      'tableInput': form.getInput('tableInput') ?? null,
    };
  }
  
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
    const customInput = await grok.functions.eval('ApiTests:CustomStringInput');
    //@ts-ignore
    form.replaceInput('stringInput', customInput);
  });

  // this test shows possible API for custom inputs of DG.InputForm
  test('replace input by custom', async () => {
    expect(inputs['stringInput'].value, 'test');
    expect(funcCall.inputs['stringInput'], 'test');
    expect(inputs['stringInput'].root.style.backgroundColor, 'aqua');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('form to funccall bind', async () => {
    inputs['stringInput'].value = 'test2';

    expect(inputs['stringInput'].value, 'test2');
    expect(funcCall.inputs['stringInput'], 'test2');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('funcall to form bind', async () => {
    funcCall.inputs['stringInput'] = 'test2';

    expect(inputs['stringInput'].value, 'test2');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('form on input change observable', async () => {
    const changedInputPropNames = [] as string[];
    const changeSub = form.onInputChanged.subscribe((ed) => {
      changedInputPropNames.push(ed.args.input.property.name);
    });

    inputs['stringInput'].value = 'test2';

    expectArray(
      changedInputPropNames,
      ['stringInput'],
    );
    changeSub.unsubscribe();
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('form to funccall bind after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;
    updateInputs();

    inputs['stringInput'].value = 'test2';

    expect(newFuncCall.inputs['stringInput'], 'test2');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('funccall to form bind after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;
    updateInputs();

    newFuncCall.inputs['stringInput'] = 'test';

    expect(inputs['stringInput'].value, 'test');
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});

  test('form on input change observable after replace', async () => {
    newFuncCall = (await grok.functions.eval('ApiTests:InputFormTest')).prepare();
    form.source = newFuncCall;
    updateInputs();

    const changedInputPropNames = [] as string[];
    form.onInputChanged.pipe(take(1)).subscribe((ed) => {
      changedInputPropNames.push(ed.args.input.property.name);
    });

    inputs['stringInput'].value = 'test2';

    expectArray(
      changedInputPropNames,
      ['stringInput'],
    );
  }, {skipReason: 'https://reddata.atlassian.net/browse/GROK-15737', owner: zakievAufar});
}, {owner: 'dkovalyov@datagrok.ai'});
