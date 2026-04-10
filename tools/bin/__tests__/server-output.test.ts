import {describe, it, expect} from 'vitest';
import {cellStr, csvCell, getKeys} from '../utils/server-output';

describe('cellStr', () => {
  it('returns empty string for null', () => {
    expect(cellStr(null)).toBe('');
  });

  it('returns empty string for undefined', () => {
    expect(cellStr(undefined)).toBe('');
  });

  it('returns the name property for objects with a name', () => {
    expect(cellStr({name: 'Alice', id: '1'})).toBe('Alice');
  });

  it('falls back to id when name is absent', () => {
    expect(cellStr({id: 'abc-123'})).toBe('abc-123');
  });

  it('returns truncated JSON for objects with neither name nor id', () => {
    const result = cellStr({foo: 'bar', baz: 42});
    expect(result).toBe('{"foo":"bar","baz":42}');
    expect(result.length).toBeLessThanOrEqual(40);
  });

  it('truncates long JSON to 40 characters', () => {
    const obj = {key: 'a'.repeat(50)};
    expect(cellStr(obj).length).toBe(40);
  });

  it('converts a number to string', () => {
    expect(cellStr(123)).toBe('123');
  });

  it('converts a boolean to string', () => {
    expect(cellStr(true)).toBe('true');
    expect(cellStr(false)).toBe('false');
  });

  it('returns a plain string as-is', () => {
    expect(cellStr('hello')).toBe('hello');
  });
});

describe('csvCell', () => {
  it('returns the string unchanged when no special chars', () => {
    expect(csvCell('hello')).toBe('hello');
  });

  it('wraps in quotes when the value contains a comma', () => {
    expect(csvCell('a,b')).toBe('"a,b"');
  });

  it('wraps and escapes double quotes', () => {
    expect(csvCell('say "hi"')).toBe('"say ""hi"""');
  });

  it('wraps in quotes when the value contains a newline', () => {
    expect(csvCell('line1\nline2')).toBe('"line1\nline2"');
  });

  it('handles an empty string', () => {
    expect(csvCell('')).toBe('');
  });
});

describe('getKeys', () => {
  it('returns an empty array for an empty input', () => {
    expect(getKeys([])).toEqual([]);
  });

  it('returns keys from a single object', () => {
    expect(getKeys([{a: 1, b: 2}])).toEqual(['a', 'b']);
  });

  it('unions keys across multiple objects', () => {
    expect(getKeys([{a: 1}, {b: 2}, {a: 3, c: 4}])).toEqual(['a', 'b', 'c']);
  });

  it('deduplicates keys', () => {
    expect(getKeys([{a: 1, b: 2}, {a: 3, b: 4}])).toEqual(['a', 'b']);
  });

  it('preserves insertion order (first appearance wins)', () => {
    expect(getKeys([{b: 1, a: 2}, {a: 3, c: 4}])).toEqual(['b', 'a', 'c']);
  });

  it('caps the result at 12 keys', () => {
    const row = Object.fromEntries(Array.from({length: 15}, (_, i) => [`k${i}`, i]));
    const keys = getKeys([row]);
    expect(keys).toHaveLength(12);
    expect(keys).toEqual(Array.from({length: 12}, (_, i) => `k${i}`));
  });
});
