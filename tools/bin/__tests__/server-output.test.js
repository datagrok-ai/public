"use strict";

var _vitest = require("vitest");
var _serverOutput = require("../utils/server-output");
(0, _vitest.describe)('cellStr', () => {
  (0, _vitest.it)('returns empty string for null', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)(null)).toBe('');
  });
  (0, _vitest.it)('returns empty string for undefined', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)(undefined)).toBe('');
  });
  (0, _vitest.it)('returns the name property for objects with a name', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)({
      name: 'Alice',
      id: '1'
    })).toBe('Alice');
  });
  (0, _vitest.it)('falls back to id when name is absent', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)({
      id: 'abc-123'
    })).toBe('abc-123');
  });
  (0, _vitest.it)('returns truncated JSON for objects with neither name nor id', () => {
    const result = (0, _serverOutput.cellStr)({
      foo: 'bar',
      baz: 42
    });
    (0, _vitest.expect)(result).toBe('{"foo":"bar","baz":42}');
    (0, _vitest.expect)(result.length).toBeLessThanOrEqual(40);
  });
  (0, _vitest.it)('truncates long JSON to 40 characters', () => {
    const obj = {
      key: 'a'.repeat(50)
    };
    (0, _vitest.expect)((0, _serverOutput.cellStr)(obj).length).toBe(40);
  });
  (0, _vitest.it)('converts a number to string', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)(123)).toBe('123');
  });
  (0, _vitest.it)('converts a boolean to string', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)(true)).toBe('true');
    (0, _vitest.expect)((0, _serverOutput.cellStr)(false)).toBe('false');
  });
  (0, _vitest.it)('returns a plain string as-is', () => {
    (0, _vitest.expect)((0, _serverOutput.cellStr)('hello')).toBe('hello');
  });
});
(0, _vitest.describe)('csvCell', () => {
  (0, _vitest.it)('returns the string unchanged when no special chars', () => {
    (0, _vitest.expect)((0, _serverOutput.csvCell)('hello')).toBe('hello');
  });
  (0, _vitest.it)('wraps in quotes when the value contains a comma', () => {
    (0, _vitest.expect)((0, _serverOutput.csvCell)('a,b')).toBe('"a,b"');
  });
  (0, _vitest.it)('wraps and escapes double quotes', () => {
    (0, _vitest.expect)((0, _serverOutput.csvCell)('say "hi"')).toBe('"say ""hi"""');
  });
  (0, _vitest.it)('wraps in quotes when the value contains a newline', () => {
    (0, _vitest.expect)((0, _serverOutput.csvCell)('line1\nline2')).toBe('"line1\nline2"');
  });
  (0, _vitest.it)('handles an empty string', () => {
    (0, _vitest.expect)((0, _serverOutput.csvCell)('')).toBe('');
  });
});
(0, _vitest.describe)('getKeys', () => {
  (0, _vitest.it)('returns an empty array for an empty input', () => {
    (0, _vitest.expect)((0, _serverOutput.getKeys)([])).toEqual([]);
  });
  (0, _vitest.it)('returns keys from a single object', () => {
    (0, _vitest.expect)((0, _serverOutput.getKeys)([{
      a: 1,
      b: 2
    }])).toEqual(['a', 'b']);
  });
  (0, _vitest.it)('unions keys across multiple objects', () => {
    (0, _vitest.expect)((0, _serverOutput.getKeys)([{
      a: 1
    }, {
      b: 2
    }, {
      a: 3,
      c: 4
    }])).toEqual(['a', 'b', 'c']);
  });
  (0, _vitest.it)('deduplicates keys', () => {
    (0, _vitest.expect)((0, _serverOutput.getKeys)([{
      a: 1,
      b: 2
    }, {
      a: 3,
      b: 4
    }])).toEqual(['a', 'b']);
  });
  (0, _vitest.it)('preserves insertion order (first appearance wins)', () => {
    (0, _vitest.expect)((0, _serverOutput.getKeys)([{
      b: 1,
      a: 2
    }, {
      a: 3,
      c: 4
    }])).toEqual(['b', 'a', 'c']);
  });
  (0, _vitest.it)('caps the result at 12 keys', () => {
    const row = Object.fromEntries(Array.from({
      length: 15
    }, (_, i) => [`k${i}`, i]));
    const keys = (0, _serverOutput.getKeys)([row]);
    (0, _vitest.expect)(keys).toHaveLength(12);
    (0, _vitest.expect)(keys).toEqual(Array.from({
      length: 12
    }, (_, i) => `k${i}`));
  });
});
(0, _vitest.describe)('printBatchOutput', () => {
  let logLines;
  let errLines;
  (0, _vitest.beforeEach)(() => {
    logLines = [];
    errLines = [];
    _vitest.vi.spyOn(console, 'log').mockImplementation((...args) => logLines.push(args.join(' ')));
    _vitest.vi.spyOn(process.stderr, 'write').mockImplementation(chunk => {
      errLines.push(String(chunk));
      return true;
    });
  });
  (0, _vitest.afterEach)(() => {
    _vitest.vi.restoreAllMocks();
  });
  (0, _vitest.it)('json: output is valid JSON equal to the full response', () => {
    const resp = {
      summary: {
        total: 1,
        succeeded: 1,
        failed: 0,
        partial: 0,
        skipped: 0
      },
      results: [{
        id: 'op0',
        action: 'users.delete',
        status: 'success'
      }]
    };
    (0, _serverOutput.printBatchOutput)(resp, 'json');
    (0, _vitest.expect)(JSON.parse(logLines[0])).toEqual(resp);
  });
  (0, _vitest.it)('table: failed operations are written to stderr as structured JSON', () => {
    const resp = {
      summary: {
        total: 1,
        succeeded: 0,
        failed: 1,
        partial: 0,
        skipped: 0
      },
      results: [{
        id: 'op0',
        action: 'users.delete',
        status: 'error',
        error: {
          error: 'Not found'
        }
      }]
    };
    (0, _serverOutput.printBatchOutput)(resp, 'table');
    (0, _vitest.expect)(errLines.length).toBeGreaterThan(0);
    const parsed = JSON.parse(errLines[0]);
    (0, _vitest.expect)(parsed[0].id).toBe('op0');
    (0, _vitest.expect)(parsed[0].error.error).toBe('Not found');
  });
});