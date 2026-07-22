/** FuncCall model, mirroring datagrok-celery-task/func_call.py.
 *  The serialized shape is produced by the Dart FuncCall.toJson() (ddt) and consumed
 *  back by datlas' `new FuncCall.fromJson(...)` (call_queue_service.dart), so toJson()
 *  must echo `func`, `options` and `aux` untouched. */

export const Type = {
  INT: 'int',
  STRING: 'string',
  BOOL: 'bool',
  FLOAT: 'double',
  BIG_INT: 'bigint',
  DATA_FRAME: 'dataframe',
  DATE_TIME: 'datetime',
  FILE: 'file',
  BLOB: 'blob',
  GRAPHICS: 'graphics',
} as const;

export const FuncCallStatus = {
  RUNNING: 'Running',
  QUEUED: 'Queued',
  COMPLETED: 'Completed',
  ERROR: 'Error',
  CANCELED: 'Canceled',
  CANCELING: 'Canceling',
} as const;

export function isCompletedStatus(status: string): boolean {
  return status === FuncCallStatus.COMPLETED || status === FuncCallStatus.ERROR ||
    status === FuncCallStatus.CANCELED;
}

export class FuncCallParam {
  name: string;
  propertyType: string;
  isInput: boolean;
  value: any;
  originalValue: any;

  constructor(name: string, propertyType: string, isInput: boolean, value: any) {
    this.name = name;
    this.propertyType = propertyType;
    this.isInput = isInput;
    this.value = value;
    this.originalValue = value;
  }

  /** The python lib (func_call.py) treats dataframe/blob/file as streamable; the Dart
   *  server (ddt func_param.dart `isStreamable`) additionally includes graphics.
   *  We match the server so `requiresPipe` agrees with datlas' decision to open a pipe. */
  get isStreamable(): boolean {
    return this.propertyType === Type.DATA_FRAME || this.propertyType === Type.BLOB ||
      this.propertyType === Type.FILE || this.propertyType === Type.GRAPHICS;
  }
}

export class FuncCall {
  static DEFAULT_BINARY_BATCH_SIZE = 2048000;

  id: string;
  func: {[key: string]: any};
  options: {[key: string]: any};
  aux: {[key: string]: any};
  params: FuncCallParam[];
  status: string;
  errorMessage: string | null = null;
  errorStackTrace: string | null = null;

  constructor(callJson: {[key: string]: any}) {
    this.id = callJson['id'];
    this.func = callJson['func'] ?? {};
    this.options = callJson['options'] ?? {};
    this.aux = callJson['aux'] ?? {};
    this.status = FuncCallStatus.RUNNING;
    const parameterValues: {[key: string]: any} = callJson['parameterValues'] ?? {};
    const funcParams: any[] = this.func['params'] ?? [];
    this.params = funcParams.map((p) => new FuncCallParam(
      p['name'], p['propertyType'], p['isInput'] ?? true, parameterValues[p['name']] ?? null));
  }

  /** Parses a celery task message body. Accepts BOTH shapes:
   *  - the Dart hybrid `{"args": [callJson], "kwargs": {}}` produced by
   *    DockerFunc.getQueueMessageBody (docker_func.dart);
   *  - celery protocol 2 `[[callJson], {}, {...embed}]` (what a real celery
   *    client would send). */
  static fromCeleryBody(body: any): FuncCall {
    const callJson = Array.isArray(body) ? body[0]?.[0] : body?.['args']?.[0];
    if (callJson == null || typeof callJson !== 'object')
      throw new Error('Unrecognized celery task body: expected {"args": [call]} or [[call], {}, {}]');
    return new FuncCall(callJson);
  }

  get funcName(): string {
    return this.func['name'] ?? '';
  }

  get useParquetTransfer(): boolean {
    return this.options['isParquet'] === true;
  }

  get userApiKey(): string | null {
    return this.aux['USER_API_KEY'] ?? null;
  }

  get apiUrl(): string | null {
    return this.aux['DATAGROK_API_URL'] ?? null;
  }

  get binaryBatchSize(): number {
    return this.aux['batchSize'] ?? FuncCall.DEFAULT_BINARY_BATCH_SIZE;
  }

  get inputParams(): FuncCallParam[] {
    return this.params.filter((p) => p.isInput);
  }

  /** Declared outputs. */
  get outputParams(): FuncCallParam[] {
    return this.params.filter((p) => !p.isInput);
  }

  get requiresPipe(): boolean {
    return this.params.some((p) => p.isStreamable);
  }

  /** Mirrors func_call.py to_json(): inputs echo the original raw values, outputs carry
   *  the serialized value (dataframe -> {id: uuid}, blob -> param name, scalars as-is). */
  toJson(): {[key: string]: any} {
    const parameterValues: {[key: string]: any} = {};
    for (const p of this.params)
      parameterValues[p.name] = p.isInput ? p.originalValue : p.value;
    return {
      'id': this.id,
      'parameterValues': parameterValues,
      'func': this.func,
      'options': this.options,
      'aux': this.aux,
      'status': this.status,
      'errorMessage': this.errorMessage,
      'errorStackTrace': this.errorStackTrace,
    };
  }
}
