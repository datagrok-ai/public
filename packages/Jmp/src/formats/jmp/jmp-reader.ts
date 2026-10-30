import {JmpFormatError} from './jmp-types';

/** Little-endian cursor over a byte buffer. Every read is bounds-checked. */
export class ByteReader {
  readonly bytes: Uint8Array;
  private readonly view: DataView;
  pos = 0;

  constructor(bytes: Uint8Array) {
    this.bytes = bytes;
    this.view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  }

  get length(): number {return this.bytes.length;}
  get eof(): boolean {return this.pos >= this.bytes.length;}

  seek(pos: number): void {
    if (pos < 0 || pos > this.bytes.length)
      throw new JmpFormatError(`Seek to ${pos} is outside the file (${this.bytes.length} bytes); truncated?`);
    this.pos = pos;
  }

  skip(n: number): void {this.seek(this.pos + n);}

  private need(n: number): void {
    if (this.pos + n > this.bytes.length)
      throw new JmpFormatError(`Read of ${n} bytes at ${this.pos} runs past end of file; file appears truncated`);
  }

  u8(): number {this.need(1); return this.view.getUint8(this.pos++);}
  i8(): number {this.need(1); return this.view.getInt8(this.pos++);}
  u16(): number {this.need(2); const v = this.view.getUint16(this.pos, true); this.pos += 2; return v;}
  i16(): number {this.need(2); const v = this.view.getInt16(this.pos, true); this.pos += 2; return v;}
  u32(): number {this.need(4); const v = this.view.getUint32(this.pos, true); this.pos += 4; return v;}
  i32(): number {this.need(4); const v = this.view.getInt32(this.pos, true); this.pos += 4; return v;}
  f64(): number {this.need(8); const v = this.view.getFloat64(this.pos, true); this.pos += 8; return v;}

  /** Reads a 64-bit integer; throws if it doesn't fit in a safe JS number. */
  i64(): number {
    this.need(8);
    const v = this.view.getBigInt64(this.pos, true);
    this.pos += 8;
    if (v > BigInt(Number.MAX_SAFE_INTEGER) || v < BigInt(Number.MIN_SAFE_INTEGER))
      throw new JmpFormatError(`64-bit value ${v} at ${this.pos - 8} is out of range`);
    return Number(v);
  }

  u64(): number {return this.i64();}

  bytesAt(n: number): Uint8Array {
    this.need(n);
    const out = this.bytes.subarray(this.pos, this.pos + n);
    this.pos += n;
    return out;
  }

  /** Reads a string prefixed by a signed length of `width` bytes (1, 2 or 4). */
  lengthPrefixedString(width: 1 | 2 | 4): string {
    const len = width === 1 ? this.i8() : width === 2 ? this.i16() : this.i32();
    if (len < 0)
      throw new JmpFormatError(`Negative string length ${len} at ${this.pos - width}`);
    return new TextDecoder('utf-8').decode(this.bytesAt(len));
  }

  /** Returns the position just after the next occurrence of `seq`, or -1 if absent. */
  indexAfter(seq: number[], from = this.pos): number {
    const b = this.bytes;
    outer: for (let i = from; i <= b.length - seq.length; i++) {
      for (let j = 0; j < seq.length; j++) {
        if (b[i + j] !== seq[j])
          continue outer;
      }
      return i + seq.length;
    }
    return -1;
  }
}
