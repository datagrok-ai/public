import {FuncCall, FuncCallParam, Type} from '../func-call';
import {marshalInput, marshalOutput, validateCall} from '../marshal';

function param(propertyType: string, value: any, isInput: boolean = true): FuncCallParam {
  return new FuncCallParam('p', propertyType, isInput, value);
}

describe('validateCall', () => {
  test('fails fast on isParquet', () => {
    const call = new FuncCall({'id': 'c', 'func': {'params': []}, 'options': {'isParquet': true}});
    expect(() => validateCall(call)).toThrow(/Parquet transfer is not supported/);
  });

  test('rejects multiple outputs (python parity message)', () => {
    const call = new FuncCall({'id': 'c', 'func': {'params': [
      {'name': 'a', 'propertyType': 'int', 'isInput': false},
      {'name': 'b', 'propertyType': 'int', 'isInput': false}]}});
    expect(() => validateCall(call)).toThrow(/Only one return parameter is allowed/);
  });

  test('accepts a single output without parquet', () => {
    const call = new FuncCall({'id': 'c', 'func': {'params': [
      {'name': 'a', 'propertyType': 'int', 'isInput': false}]}, 'options': {}});
    expect(() => validateCall(call)).not.toThrow();
  });
});

describe('marshalInput', () => {
  test('int: number and numeric string', () => {
    expect(marshalInput(param(Type.INT, 7))).toBe(7);
    expect(marshalInput(param(Type.INT, '7'))).toBe(7);
    expect(() => marshalInput(param(Type.INT, 'abc'))).toThrow(/Incorrect input type/);
    expect(() => marshalInput(param(Type.INT, 1.5))).toThrow(/Incorrect input type/);
  });

  test('double', () => {
    expect(marshalInput(param(Type.FLOAT, 1.5))).toBe(1.5);
    expect(marshalInput(param(Type.FLOAT, '1.5'))).toBe(1.5);
    expect(() => marshalInput(param(Type.FLOAT, 'abc'))).toThrow(/Incorrect input type/);
  });

  test('bool', () => {
    expect(marshalInput(param(Type.BOOL, true))).toBe(true);
    expect(marshalInput(param(Type.BOOL, 'false'))).toBe(false);
    expect(() => marshalInput(param(Type.BOOL, 'yes'))).toThrow(/Incorrect input type/);
  });

  test('string', () => {
    expect(marshalInput(param(Type.STRING, 'hi'))).toBe('hi');
    expect(marshalInput(param(Type.STRING, 5))).toBe('5');
  });

  test('bigint via BigInt(string)', () => {
    expect(marshalInput(param(Type.BIG_INT, '9007199254740993'))).toBe(BigInt('9007199254740993'));
    expect(() => marshalInput(param(Type.BIG_INT, 'nope'))).toThrow(/bigint/);
  });

  test('datetime from ISO string', () => {
    const value = marshalInput(param(Type.DATE_TIME, '2026-01-02T03:04:05.000Z'));
    expect(value).toBeInstanceOf(Date);
    expect((value as Date).toISOString()).toBe('2026-01-02T03:04:05.000Z');
    expect(() => marshalInput(param(Type.DATE_TIME, 'not-a-date'))).toThrow(/ISO datetime/);
  });

  test('dataframe: UTF-8 CSV bytes through DG.DataFrame.fromCsv', () => {
    const fromCsv = jest.fn().mockReturnValue({'fake': 'df'});
    const dg = {DataFrame: {fromCsv: fromCsv}};
    const p = param(Type.DATA_FRAME, new Uint8Array(Buffer.from('a,b\n1,2\n', 'utf8')));
    expect(marshalInput(p, dg)).toEqual({'fake': 'df'});
    expect(fromCsv).toHaveBeenCalledWith('a,b\n1,2\n');
    expect(() => marshalInput(param(Type.DATA_FRAME, 'not-bytes'), dg)).toThrow(/Incorrect input type/);
  });

  test('blob passes bytes through', () => {
    const bytes = new Uint8Array([1, 2, 3]);
    expect(marshalInput(param(Type.BLOB, bytes))).toBe(bytes);
    expect(() => marshalInput(param(Type.BLOB, 'no'))).toThrow(/Incorrect input type/);
  });

  test('null stays null', () => {
    expect(marshalInput(param(Type.INT, null))).toBeNull();
  });
});

describe('marshalOutput', () => {
  test('dataframe: csv bytes + value {id} + tags', () => {
    const p = param(Type.DATA_FRAME, null, false);
    const df = {toCsv: () => 'a,b\n1,2\n'};
    const out = marshalOutput(p, df);
    expect(Buffer.from(out.bytes!).toString('utf8')).toBe('a,b\n1,2\n');
    expect(p.value).toEqual({'id': out.tags!['.id']});
    expect(out.tags!['.type']).toBe('csv');
    expect(() => marshalOutput(param(Type.DATA_FRAME, null, false), 'not-a-df')).toThrow(/Incorrect return type/);
  });

  test('blob: bytes + value = param name + tags', () => {
    const p = param(Type.BLOB, null, false);
    const out = marshalOutput(p, new Uint8Array([9, 8]));
    expect(Array.from(out.bytes!)).toEqual([9, 8]);
    expect(p.value).toBe('p');
    expect(out.tags).toEqual({'.id': 'p', '.type': 'blob'});
    expect(() => marshalOutput(param(Type.BLOB, null, false), 'str')).toThrow(/Incorrect return type/);
  });

  test('scalars coerced, NaN throws', () => {
    const intParam = param(Type.INT, null, false);
    marshalOutput(intParam, 3.7);
    expect(intParam.value).toBe(3); // python int() truncation parity
    expect(() => marshalOutput(param(Type.INT, null, false), 'abc')).toThrow(/Incorrect return type/);
    const floatParam = param(Type.FLOAT, null, false);
    marshalOutput(floatParam, '2.5');
    expect(floatParam.value).toBe(2.5);
    expect(() => marshalOutput(param(Type.FLOAT, null, false), 'abc')).toThrow(/Incorrect return type/);
    const boolParam = param(Type.BOOL, null, false);
    marshalOutput(boolParam, 1);
    expect(boolParam.value).toBe(true);
    const strParam = param(Type.STRING, null, false);
    marshalOutput(strParam, 42);
    expect(strParam.value).toBe('42');
    const bigParam = param(Type.BIG_INT, null, false);
    marshalOutput(bigParam, BigInt('9007199254740993'));
    expect(bigParam.value).toBe('9007199254740993');
  });

  test('datetime output requires Date/dayjs', () => {
    const p = param(Type.DATE_TIME, null, false);
    marshalOutput(p, new Date('2026-01-02T03:04:05.000Z'));
    expect(p.value).toBe('2026-01-02T03:04:05.000Z');
    expect(() => marshalOutput(param(Type.DATE_TIME, null, false), '2026-01-02')).toThrow(/Incorrect return type/);
  });

  test('null return keeps value null, no streaming', () => {
    const p = param(Type.DATA_FRAME, null, false);
    const out = marshalOutput(p, null);
    expect(p.value).toBeNull();
    expect(out.bytes).toBeNull();
    expect(out.tags).toBeNull();
  });

  test('unsupported return type throws', () => {
    expect(() => marshalOutput(param(Type.GRAPHICS, null, false), 'svg')).toThrow(/Unsupported return type/);
  });
});
