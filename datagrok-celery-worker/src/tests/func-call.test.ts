import {FuncCall, FuncCallStatus} from '../func-call';

function sampleCallJson(): any {
  return {
    'id': 'call-1',
    'func': {
      'name': 'myFunc',
      'params': [
        {'name': 'x', 'propertyType': 'int', 'isInput': true},
        {'name': 'df', 'propertyType': 'dataframe', 'isInput': true},
        {'name': 'result', 'propertyType': 'dataframe', 'isInput': false},
      ],
    },
    'parameterValues': {'x': 42},
    'options': {'isParquet': false},
    'aux': {'USER_API_KEY': 'secret', 'DATAGROK_API_URL': 'http://localhost:8082/api', 'batchSize': 1000},
  };
}

describe('FuncCall', () => {
  test('parses the Dart hybrid body {"args":[call],"kwargs":{}}', () => {
    const call = FuncCall.fromCeleryBody({'args': [sampleCallJson()], 'kwargs': {}});
    expect(call.id).toBe('call-1');
    expect(call.funcName).toBe('myFunc');
    expect(call.params.length).toBe(3);
    expect(call.params[0].value).toBe(42);
    expect(call.status).toBe(FuncCallStatus.RUNNING);
  });

  test('parses the celery proto2 array body [[call],{},{}]', () => {
    const call = FuncCall.fromCeleryBody([[sampleCallJson()], {}, {}]);
    expect(call.id).toBe('call-1');
    expect(call.inputParams.map((p) => p.name)).toEqual(['x', 'df']);
    expect(call.outputParams.map((p) => p.name)).toEqual(['result']);
  });

  test('throws on an unrecognized body', () => {
    expect(() => FuncCall.fromCeleryBody({'kwargs': {}})).toThrow(/Unrecognized celery task body/);
    expect(() => FuncCall.fromCeleryBody(null)).toThrow(/Unrecognized celery task body/);
  });

  test('streamable detection: dataframe/blob/file/graphics stream, scalars do not', () => {
    const call = new FuncCall({
      'id': 'c',
      'func': {'name': 'f', 'params': [
        {'name': 'a', 'propertyType': 'int', 'isInput': true},
        {'name': 'b', 'propertyType': 'blob', 'isInput': true},
        {'name': 'c', 'propertyType': 'file', 'isInput': true},
        {'name': 'g', 'propertyType': 'graphics', 'isInput': false},
      ]},
    });
    expect(call.params.map((p) => p.isStreamable)).toEqual([false, true, true, true]);
    expect(call.requiresPipe).toBe(true);
    const scalarCall = new FuncCall({'id': 'c2', 'func': {'name': 'f', 'params': [
      {'name': 'a', 'propertyType': 'string', 'isInput': true}]}});
    expect(scalarCall.requiresPipe).toBe(false);
  });

  test('aux getters', () => {
    const call = new FuncCall(sampleCallJson());
    expect(call.userApiKey).toBe('secret');
    expect(call.apiUrl).toBe('http://localhost:8082/api');
    expect(call.binaryBatchSize).toBe(1000);
    expect(new FuncCall({'id': 'x'}).binaryBatchSize).toBe(FuncCall.DEFAULT_BINARY_BATCH_SIZE);
  });

  test('toJson round-trip: aux/func/options echoed, inputs keep original values, df output replaced by {id}', () => {
    const source = sampleCallJson();
    const call = new FuncCall(source);
    call.params[1].value = new Uint8Array([1, 2, 3]); // streamed + marshaled input
    call.params[2].value = {'id': 'some-uuid'};       // marshaled df output
    call.status = FuncCallStatus.COMPLETED;
    const json = call.toJson();
    expect(json['id']).toBe('call-1');
    expect(json['func']).toEqual(source['func']);
    expect(json['options']).toEqual(source['options']);
    expect(json['aux']).toEqual(source['aux']);
    expect(json['status']).toBe('Completed');
    expect(json['errorMessage']).toBeNull();
    expect(json['errorStackTrace']).toBeNull();
    // inputs echo the ORIGINAL raw values, not the marshaled ones
    expect(json['parameterValues']['x']).toBe(42);
    expect(json['parameterValues']['df']).toBeNull();
    expect(json['parameterValues']['result']).toEqual({'id': 'some-uuid'});
  });
});
