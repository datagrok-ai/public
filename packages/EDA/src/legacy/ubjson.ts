const OBJECT = 0x7b;
const OBJECT_END = 0x7d;
const ARRAY = 0x5b;
const ARRAY_END = 0x5d;
const TYPE = 0x24;
const COUNT = 0x23;
const INT64 = 0x4c;

/** JSON text of a UBJSON value in the flavour XGBoost writes (`L`-typed lengths, typed arrays `[$t#L n`, counted
 * arrays `[#L n`); floats keep a fraction or an exponent, so XGBoost's text reader reads them as numbers, integers
 * stay bare */
export function ubjsonToJson(bytes: Uint8Array): string {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const decoder = new TextDecoder();
  let pos = 0;

  const unsupported = (at: number): Error => new Error(`UBJSON: unsupported marker ${bytes[at]} at ${at}`);

  const read = <T>(size: number, get: (at: number) => T): T => {
    const v = get(pos);
    pos += size;
    return v;
  };

  const length = (): number => {
    if (bytes[pos] !== INT64)
      throw unsupported(pos);
    pos++;
    const n = Number(read(8, (at) => view.getBigInt64(at)));
    if (n < 0 || n > bytes.length - pos)
      throw new Error(`UBJSON: length ${n} out of range at ${pos - 9}`);
    return n;
  };

  const text = (): string => {
    const n = length();
    return read(n, (at) => decoder.decode(bytes.subarray(at, at + n)));
  };

  const float = (v: number): string => {
    if (!Number.isFinite(v))
      throw new Error(`UBJSON: non-finite float ${v} at ${pos}`);
    const s = String(v);
    return Number.isInteger(v) && Math.abs(v) < 1e21 ? `${s}.0` : s;
  };

  const scalar = (marker: number, at: number): string => {
    switch (marker) {
    case 0x5a: return 'null';
    case 0x54: return 'true';
    case 0x46: return 'false';
    case 0x64: return float(read(4, (i) => view.getFloat32(i)));
    case 0x44: return float(read(8, (i) => view.getFloat64(i)));
    case 0x53: return JSON.stringify(text());
    case 0x69: return String(read(1, (i) => view.getInt8(i)));
    case 0x55: return String(read(1, (i) => view.getUint8(i)));
    case 0x49: return String(read(2, (i) => view.getInt16(i)));
    case 0x6c: return String(read(4, (i) => view.getInt32(i)));
    case INT64: return String(read(8, (i) => view.getBigInt64(i)));
    case 0x43: return JSON.stringify(String.fromCharCode(read(1, (i) => view.getUint8(i))));
    default: throw unsupported(at);
    }
  };

  /** `n` items, or items up to the closing marker of the container when `n` is undefined, comma-separated */
  const seq = (item: () => string, n?: number): string => {
    const items: string[] = [];
    if (n === undefined) {
      while (bytes[pos] !== OBJECT_END && bytes[pos] !== ARRAY_END)
        items.push(item());
      pos++;
    } else {
      for (let i = 0; i < n; ++i)
        items.push(item());
    }
    return items.join(',');
  };

  const value = (): string => {
    const at = pos;
    const marker = bytes[pos++];
    if (marker === OBJECT) {
      if (bytes[pos] === TYPE || bytes[pos] === COUNT)
        throw unsupported(pos);
      return `{${seq(() => `${JSON.stringify(text())}:${value()}`)}}`;
    }
    if (marker !== ARRAY)
      return scalar(marker, at);
    if (bytes[pos] === TYPE) {
      const type = bytes[pos + 1];
      if (bytes[pos + 2] !== COUNT)
        throw unsupported(pos + 2);
      pos += 3;
      return `[${seq(() => scalar(type, pos), length())}]`;
    }
    if (bytes[pos] === COUNT) {
      pos++;
      return `[${seq(value, length())}]`;
    }
    return `[${seq(value)}]`;
  };

  return value();
}
