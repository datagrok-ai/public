import * as DG from 'datagrok-api/dg';

// Support for json serialization/deserialization of objects that
// include DataFrames. Not efficient, intended for writing tests.

const customTypeKey = '_DG_CUSTOM_SERIALIZED_TOKEN_';

export function serialize(obj: any, space = 2) {
  return JSON.stringify(
    obj,
    (_key, value) => {
      if (value instanceof DG.DataFrame) {
        return {
          [customTypeKey]: 'DataFrame',
          value: Array.from(value.toByteArray())
        };
      }
      if (value instanceof ArrayBuffer) {
        return {
          [customTypeKey]: 'ArrayBuffer',
          value: Array.from(new Uint8Array(value))
        };
      }
      if (value instanceof Map) {
        return {
          [customTypeKey]: 'Map',
          value: Array.from(value)
        };
      }
      if (value instanceof Set) {
        return {
          [customTypeKey]: 'Set',
          value: Array.from(value)
        };
      }
      return value;
    },
    space
  );
}

export function deserialize(obj: string) {
  return JSON.parse(
    obj,
    transform
  );
}

export function applyTransformations(obj: any) {
  return deserialize(JSON.stringify(obj));
}

// just apply data transformations if needed
export function transform(_key: string, value: any) {
  if (value && value[customTypeKey] && value.value) {
    switch (value[customTypeKey]) {
    case 'DataFrame':
      return DG.DataFrame.fromByteArray(new Uint8Array(value.value));
    case 'ArrayBuffer':
      return new Uint8Array(value.value).buffer;
    case 'Map':
      return new Map(value.value);
    case 'Set':
      return new Set(value.value);
    }
  }
  return value;
}
