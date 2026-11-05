import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {irisMeasurements, mulberry32, seededRegression} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const TIMEOUT = 300000;
const DATA_INPUTS = ['df', 'predictColumn', 'seed'];
const REGRESSORS = ['Linear Regression', 'PLS Regression'];
const CLASSIFIERS = ['Softmax', 'Naive Bayes', 'LDA'];
const REGRESSION_ROWS = 30;
const MAX_REGRESSION_EPOCHS = 64;
const POINTS = 5;

type Value = number | string | boolean;
type Dimension = {type: string, low?: number, high?: number, log?: boolean, step?: number, choices?: string[],
  simpler?: string};
type Space = {[name: string]: Dimension};
type Data = {task: string, df: DG.DataFrame, predictColumn: DG.Column};

function trainFunctions(): DG.Func[] {
  return DG.Func.find({package: 'Eda', meta: {mlrole: 'train'}});
}

function spaceOf(f: DG.Func): Space | undefined {
  const text = f.options['mlspace'];
  return text === undefined || text === null ? undefined : JSON.parse(text);
}

function defaultsOf(f: DG.Func): {[name: string]: Value} {
  const values: {[name: string]: Value} = {};
  for (const p of f.inputs.filter((p) => !['df', 'predictColumn'].includes(p.name))) {
    const text = p.initialValue;
    if (typeof text === 'string' && text !== '') {
      values[p.name] = p.propertyType === DG.TYPE.STRING ? text.replace(/^(["'])(.*)\1$/, '$2') :
        p.propertyType === DG.TYPE.BOOL ? text.toLowerCase() === 'true' : Number(text);
    }
  }
  return values;
}

function sample(dim: Dimension, rng: () => number): Value {
  const u = rng();
  if (dim.type === 'categorical')
    return dim.choices![Math.floor(u * dim.choices!.length)];
  const low = dim.low!;
  const high = dim.high!;
  if (dim.type === 'int') {
    const v = dim.log ? Math.round(Math.exp(Math.log(low) + u * (Math.log(high) - Math.log(low)))) :
      low + Math.floor(u * (high - low + 1));
    return dim.step === undefined ? v : Math.round((v - low) / dim.step) * dim.step + low;
  }
  return dim.log ? Math.exp(Math.log(low) + u * (Math.log(high) - Math.log(low))) : low + u * (high - low);
}

async function datasets(): Promise<Data[]> {
  const iris = await grok.dapi.files.readCsv(IRIS);
  const regression = seededRegression(REGRESSION_ROWS, 15, 1);
  const y = regression.col('y')!;
  regression.columns.remove('y');
  return [{task: 'classification', df: irisMeasurements(iris), predictColumn: iris.col('Species')!},
    {task: 'regression', df: regression, predictColumn: y}];
}

function dataOf(f: DG.Func, all: Data[]): Data[] {
  const method = f.options['mlname'];
  return all.filter((d) => d.task === 'regression' ? !CLASSIFIERS.includes(method) : !REGRESSORS.includes(method));
}

async function fitFailure(f: DG.Func, data: Data, values: {[name: string]: Value}, label: string):
  Promise<string | null> {
  try {
    const model = await f.apply({...values, df: data.df, predictColumn: data.predictColumn});
    return model instanceof Uint8Array && model.length > 0 ? null : `${label}: no model`;
  } catch (e) {
    return `${label}: ${e instanceof Error ? e.message : e}`;
  }
}

category('Crux search spaces', () => {
  test('Every method declares mlparallel and the spaces are well formed', async () => {
    const functions = trainFunctions();
    expect(functions.length, 16);
    const problems: string[] = [];
    let spaces = 0;
    for (const f of functions) {
      const method = f.options['mlname'];
      if (f.options['mlparallel'] !== 'true')
        problems.push(`${method}: no mlparallel`);
      const space = spaceOf(f);
      if (space === undefined)
        continue;
      ++spaces;
      for (const [name, dim] of Object.entries(space)) {
        const at = `${method} ${name}`;
        const input = f.inputs.find((p) => p.name === name);
        if (input === undefined || DATA_INPUTS.includes(name)) {
          problems.push(`${at}: not a searchable input`);
          continue;
        }
        if (!['int', 'float', 'categorical'].includes(dim.type))
          problems.push(`${at}: type ${dim.type}`);
        if (dim.simpler !== undefined && !(['low', 'high'].includes(dim.simpler) ||
          (dim.simpler === 'first' && dim.type === 'categorical')))
          problems.push(`${at}: simpler ${dim.simpler}`);
        if (dim.type === 'categorical') {
          const choices = input.choices ?? [];
          if (!Array.isArray(dim.choices) || dim.choices.length === 0)
            problems.push(`${at}: no choices`);
          else if (choices.length > 0 && dim.choices.some((c) => !choices.includes(c)))
            problems.push(`${at}: choices ${dim.choices} not in ${choices}`);
          continue;
        }
        const low = dim.low!;
        const high = dim.high!;
        if (!(low < high))
          problems.push(`${at}: low ${low} >= high ${high}`);
        if (dim.log && !(dim.type === 'int' ? low >= 1 : low > 0))
          problems.push(`${at}: log with low ${low}`);
        if (typeof input.min === 'number' && low < input.min)
          problems.push(`${at}: low ${low} below the input min ${input.min}`);
        if (typeof input.max === 'number' && high > input.max)
          problems.push(`${at}: high ${high} above the input max ${input.max}`);
      }
    }
    expect(problems.join('\n'), '');
    expect(spaces, 15);
  });

  test('Every mlenabled rule names real inputs', async () => {
    const problems: string[] = [];
    let rules = 0;
    const isObject = (v: unknown): v is {[key: string]: unknown} =>
      typeof v === 'object' && v !== null && !Array.isArray(v);
    for (const f of trainFunctions()) {
      const text = f.options['mlenabled'];
      if (text === undefined || text === null)
        continue;
      ++rules;
      const method = f.options['mlname'];
      const inputOf = (name: string) => ['df', 'predictColumn'].includes(name) ? undefined :
        f.inputs.find((p) => p.name === name);
      let parsed: unknown;
      try {
        parsed = JSON.parse(text);
      } catch {
        problems.push(`${method}: mlenabled is not JSON`);
        continue;
      }
      if (!isObject(parsed)) {
        problems.push(`${method}: mlenabled is not an object`);
        continue;
      }
      for (const [name, rule] of Object.entries(parsed)) {
        const at = `${method} ${name}`;
        if (inputOf(name) === undefined)
          problems.push(`${at}: not an input`);
        const conditions = Array.isArray(rule) ? rule : [rule];
        if (conditions.length === 0)
          problems.push(`${at}: no conditions`);
        for (const condition of conditions) {
          if (!isObject(condition) || Object.keys(condition).length === 0) {
            problems.push(`${at}: a condition is not a non-empty object`);
            continue;
          }
          for (const [other, clause] of Object.entries(condition)) {
            if (other === 'task') {
              if (!Array.isArray(clause) || clause.length === 0 ||
                clause.some((t) => !['classification', 'regression'].includes(t)))
                problems.push(`${at}: task ${JSON.stringify(clause)}`);
              continue;
            }
            const input = inputOf(other);
            if (input === undefined || other === name) {
              problems.push(`${at}: clause on ${other}, not another input`);
              continue;
            }
            if (Array.isArray(clause)) {
              const choices = input.choices ?? [];
              if (clause.length === 0 || (choices.length > 0 && clause.some((c) => !choices.includes(c))))
                problems.push(`${at}: values ${JSON.stringify(clause)} not in ${choices}`);
            } else if (!isObject(clause) || Object.keys(clause).length === 0 ||
              Object.entries(clause).some(([op, v]) => !['gt', 'lt'].includes(op) || typeof v !== 'number'))
              problems.push(`${at}: clause ${JSON.stringify(clause)}`);
          }
        }
      }
    }
    expect(problems.join('\n'), '');
    expect(rules, 3);
  });

  test('The corners of every space fit', async () => {
    const all = await datasets();
    const failures: string[] = [];
    let corners = 0;
    for (const f of trainFunctions()) {
      const space = spaceOf(f);
      if (space === undefined)
        continue;
      const defaults = defaultsOf(f);
      const method = f.options['mlname'];
      for (const data of dataOf(f, all)) {
        for (const [name, dim] of Object.entries(space)) {
          const ends: Value[] = dim.type === 'categorical' ? [dim.choices![0], dim.choices![dim.choices!.length - 1]] :
            [dim.low!, dim.high!];
          for (const value of ends) {
            const tooSlow = data.task === 'regression' && name === 'epochs' &&
              (value as number) > MAX_REGRESSION_EPOCHS;
            const tooManyNeighbors = data.task === 'regression' && name === 'neighbors' &&
              (value as number) > REGRESSION_ROWS;
            if (tooSlow || tooManyNeighbors)
              continue;
            ++corners;
            const failure = await fitFailure(f, data, {...defaults, [name]: value},
              `${method} (${data.task}) ${name}=${value}`);
            if (failure !== null)
              failures.push(failure);
          }
        }
      }
    }
    console.log(`Crux search spaces: ${corners} corners fitted`);
    expect(failures.join('\n'), '');
  }, {timeout: TIMEOUT});

  test('Five random points per space fit', async () => {
    const all = await datasets();
    const rng = mulberry32(4);
    const failures: string[] = [];
    for (const f of trainFunctions()) {
      const space = spaceOf(f);
      if (space === undefined)
        continue;
      const defaults = defaultsOf(f);
      const data = dataOf(f, all)[0];
      for (let k = 0; k < POINTS; ++k) {
        const point: {[name: string]: Value} = {};
        for (const [name, dim] of Object.entries(space))
          point[name] = sample(dim, rng);
        const failure = await fitFailure(f, data, {...defaults, ...point},
          `${f.options['mlname']} (${data.task}) ${JSON.stringify(point)}`);
        if (failure !== null)
          failures.push(failure);
      }
    }
    expect(failures.join('\n'), '');
  }, {timeout: TIMEOUT});
});
