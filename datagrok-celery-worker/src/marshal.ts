/** Input/output value marshaling, mirroring datagrok-celery-task/utils.py
 *  (InputValueProcessor/ReturnValueProcessor) with DG.DataFrame instead of pandas. */
import {randomUUID} from 'node:crypto';

import {FuncCall, FuncCallParam, Type} from './func-call';

/** Fail fast on unsupported call shapes (python DatagrokTask.__call__ parity). */
export function validateCall(call: FuncCall): void {
  // The server sets options.isParquet = false for js-lang funcs (docker_func.dart
  // executeCall); a parquet payload means a misconfigured deployment.
  if (call.useParquetTransfer)
    throw new Error('Parquet transfer is not supported by the JS worker: the server must send dataframes as d42 or CSV (call.options.isParquet must be false)');
  if (call.outputParams.length > 1)
    throw new Error('Only one return parameter is allowed, but more are present in the Datagrok function declaration.');
}

function inputError(param: FuncCallParam, expected: string): Error {
  return new Error(`Incorrect input type for ${param.name}. Expected ${expected}, got ${typeof param.value}`);
}

/** Coerces param.value (raw JSON value, or Uint8Array for streamed params) into the
 *  runtime value passed to the package function. [dg] is the datagrok-api DG namespace
 *  (required for dataframe params only, so scalar marshaling is testable standalone). */
export function marshalInput(param: FuncCallParam, dg?: any): any {
  const value = param.value;
  if (value == null)
    return null;
  switch (param.propertyType) {
    case Type.INT: {
      const n = typeof value === 'number' ? value : typeof value === 'string' ? Number(value) : NaN;
      if (!Number.isInteger(n))
        throw inputError(param, 'int');
      param.value = n;
      break;
    }
    case Type.FLOAT: {
      const n = typeof value === 'number' ? value : typeof value === 'string' ? Number(value) : NaN;
      if (Number.isNaN(n))
        throw inputError(param, 'double');
      param.value = n;
      break;
    }
    case Type.BOOL: {
      if (typeof value === 'boolean')
        break;
      if (value === 'true' || value === 'false')
        param.value = value === 'true';
      else
        throw inputError(param, 'bool');
      break;
    }
    case Type.STRING:
      param.value = String(value);
      break;
    case Type.BIG_INT: {
      try {
        param.value = BigInt(typeof value === 'string' ? value : String(value));
      }
      catch (_) {
        throw new Error(`Incorrect input type for ${param.name}. Expected a string-representable bigint.`);
      }
      break;
    }
    case Type.DATE_TIME: {
      // ISO string passthrough (the server serializes datetimes as ISO-8601)
      if (value instanceof Date)
        break;
      if (typeof value !== 'string')
        throw inputError(param, 'ISO datetime string');
      const date = new Date(value);
      if (Number.isNaN(date.getTime()))
        throw new Error(`Incorrect input value for ${param.name}: not a valid ISO datetime: ${value}`);
      param.value = date;
      break;
    }
    case Type.DATA_FRAME: {
      if (!(value instanceof Uint8Array))
        throw inputError(param, 'bytes (streamed dataframe)');
      if (dg?.DataFrame == null)
        throw new Error('DG runtime is not initialized: cannot parse a dataframe input');
      // the pipe '.type' tag picks the decoder: 'dataframe' is native d42 (what
      // server_action.dart sends to js-lang funcs), 'csv' comes from older datlas versions
      if (param.receivedType === Type.DATA_FRAME)
        param.value = dg.DataFrame.fromByteArray(value);
      else if (param.receivedType === 'csv')
        param.value = dg.DataFrame.fromCsv(Buffer.from(value.buffer, value.byteOffset, value.byteLength).toString('utf8'));
      else
        throw new Error(`Unsupported dataframe transfer type '${param.receivedType}' for ${param.name}: expected 'dataframe' (d42 binary) or 'csv'`);
      break;
    }
    case Type.FILE:
    case Type.BLOB: {
      if (!(value instanceof Uint8Array))
        throw inputError(param, 'bytes');
      break;
    }
    default:
      break; // unknown types pass through as-is (python parity: no processor registered)
  }
  return param.value;
}

export interface MarshaledOutput {
  /** Bytes to stream over the pipe for dataframe/blob outputs, null otherwise. */
  bytes: Uint8Array | null;
  /** grok_pipe SENDING tags ('.id'/'.type') for streamable outputs, null otherwise. */
  tags: {[key: string]: string} | null;
}

function returnError(param: FuncCallParam, expected: string, value: any): Error {
  return new Error(`Incorrect return type for ${param.name}. Expected ${expected}, got ${typeof value}`);
}

/** Sets param.value to the serialized value that goes into the result CALL json
 *  (mirrors ReturnValueProcessor + _send_param_grok_pipe: dataframe -> {id: uuid} with
 *  tags {'.id': id, '.type': 'dataframe'} and d42 bytes, blob -> param name with
 *  {'.id': name, '.type': 'blob'}) and returns the bytes to stream for streamable outputs. */
export function marshalOutput(param: FuncCallParam, value: any, dg?: any): MarshaledOutput {
  if (value == null) {
    param.value = null;
    return {bytes: null, tags: null};
  }
  switch (param.propertyType) {
    case Type.DATA_FRAME: {
      if (typeof value?.toByteArray !== 'function')
        throw returnError(param, 'DG.DataFrame', value);
      const id = randomUUID();
      param.value = {'id': id};
      // native d42 — the Dart receiver (sockets.dart GrokSocketDecoder) decodes any
      // non-csv dataframe with DataFrame.fromByteArray
      return {bytes: value.toByteArray(), tags: {'.id': id, '.type': Type.DATA_FRAME}};
    }
    // python's ReturnValueProcessor supports blob only; the Dart server handles FILE
    // identically to BLOB (server_action.dart getParamsFromRemoteCall), so file outputs
    // are accepted here as a superset.
    case Type.FILE:
    case Type.BLOB: {
      if (!(value instanceof Uint8Array))
        throw returnError(param, 'Uint8Array', value);
      param.value = param.name;
      return {bytes: value instanceof Buffer ? new Uint8Array(value) : value, tags: {'.id': param.name, '.type': 'blob'}};
    }
    case Type.INT: {
      const n = Number(value);
      if (Number.isNaN(n))
        throw returnError(param, 'int', value);
      param.value = Math.trunc(n);
      break;
    }
    case Type.FLOAT: {
      const n = Number(value);
      if (Number.isNaN(n))
        throw returnError(param, 'double', value);
      param.value = n;
      break;
    }
    case Type.BOOL:
      param.value = Boolean(value);
      break;
    case Type.STRING:
      param.value = String(value);
      break;
    case Type.BIG_INT:
      param.value = String(value); // python set_bigint: str(val)
      break;
    case Type.DATE_TIME: {
      if (typeof value?.toISOString !== 'function') // Date or dayjs
        throw returnError(param, 'Date/dayjs', value);
      param.value = value.toISOString();
      break;
    }
    default:
      throw new Error(`Unsupported return type: ${param.propertyType}`);
  }
  return {bytes: null, tags: null};
}
