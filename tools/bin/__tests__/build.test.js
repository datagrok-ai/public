"use strict";

var _vitest = require("vitest");
var _build = require("../commands/build");
(0, _vitest.describe)('getNestedValue', () => {
  (0, _vitest.it)('returns value for a simple key', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      name: 'Chem'
    }, 'name')).toBe('Chem');
  });
  (0, _vitest.it)('returns value for a nested path', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      a: {
        b: {
          c: 42
        }
      }
    }, 'a.b.c')).toBe(42);
  });
  (0, _vitest.it)('returns undefined for a missing key', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      name: 'Chem'
    }, 'version')).toBeUndefined();
  });
  (0, _vitest.it)('returns undefined when a mid-path segment is null', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      a: null
    }, 'a.b')).toBeUndefined();
  });
  (0, _vitest.it)('returns undefined when a mid-path segment is missing', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      a: {}
    }, 'a.b.c')).toBeUndefined();
  });
  (0, _vitest.it)('returns undefined for an empty path (splits to empty string key)', () => {
    (0, _vitest.expect)((0, _build.getNestedValue)({
      x: 1
    }, '')).toBeUndefined();
  });
});
const pkg = overrides => ({
  dir: '/tmp/pkg',
  name: overrides.name ?? 'test-pkg',
  friendlyName: overrides.friendlyName ?? overrides.name ?? 'Test Pkg',
  version: overrides.version ?? '1.0.0',
  packageJson: overrides
});
(0, _vitest.describe)('applyFilter', () => {
  const packages = [pkg({
    name: 'Chem',
    version: '1.5.0',
    category: 'Cheminformatics'
  }), pkg({
    name: 'Bio',
    version: '2.0.0',
    category: 'Bioinformatics'
  }), pkg({
    name: 'PowerGrid',
    version: '1.5.0',
    category: 'Viewers'
  })];
  (0, _vitest.it)('returns all packages when filter matches all', () => {
    (0, _vitest.expect)((0, _build.applyFilter)(packages, 'name:.')).toHaveLength(3);
  });
  (0, _vitest.it)('filters by exact name match', () => {
    const result = (0, _build.applyFilter)(packages, 'name:^Chem$');
    (0, _vitest.expect)(result).toHaveLength(1);
    (0, _vitest.expect)(result[0].name).toBe('Chem');
  });
  (0, _vitest.it)('filters by partial name (regex substring)', () => {
    const result = (0, _build.applyFilter)(packages, 'name:Bio');
    (0, _vitest.expect)(result).toHaveLength(1);
    (0, _vitest.expect)(result[0].name).toBe('Bio');
  });
  (0, _vitest.it)('returns empty array when nothing matches', () => {
    (0, _vitest.expect)((0, _build.applyFilter)(packages, 'name:NOMATCH')).toHaveLength(0);
  });
  (0, _vitest.it)('filters by version', () => {
    const result = (0, _build.applyFilter)(packages, 'version:^1\\.5');
    (0, _vitest.expect)(result).toHaveLength(2);
    (0, _vitest.expect)(result.map(p => p.name)).toEqual(_vitest.expect.arrayContaining(['Chem', 'PowerGrid']));
  });
  (0, _vitest.it)('applies && conjunction (both conditions must match)', () => {
    const result = (0, _build.applyFilter)(packages, 'name:Chem && version:1\\.5');
    (0, _vitest.expect)(result).toHaveLength(1);
    (0, _vitest.expect)(result[0].name).toBe('Chem');
  });
  (0, _vitest.it)('returns empty when one part of && conjunction fails', () => {
    (0, _vitest.expect)((0, _build.applyFilter)(packages, 'name:Chem && version:^2')).toHaveLength(0);
  });
  (0, _vitest.it)('filters by nested field', () => {
    const withNested = [pkg({
      name: 'A',
      datagrok: {
        apiVersion: '1.0'
      }
    }), pkg({
      name: 'B',
      datagrok: {
        apiVersion: '2.0'
      }
    })];
    const result = (0, _build.applyFilter)(withNested, 'datagrok.apiVersion:^1');
    (0, _vitest.expect)(result).toHaveLength(1);
    (0, _vitest.expect)(result[0].name).toBe('A');
  });
  (0, _vitest.it)('returns empty when field does not exist', () => {
    (0, _vitest.expect)((0, _build.applyFilter)(packages, 'nonexistent:anything')).toHaveLength(0);
  });
  (0, _vitest.it)('treats filter with no colon as field name with match-all pattern', () => {
    // No colon → field = whole string, pattern = /./  (matches any value)
    // The function returns packages where the field exists and is non-empty
    const result = (0, _build.applyFilter)(packages, 'name');
    (0, _vitest.expect)(result).toHaveLength(3);
  });
});