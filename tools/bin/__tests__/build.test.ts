import {describe, it, expect} from 'vitest';
import {getNestedValue, applyFilter} from '../commands/build';

describe('getNestedValue', () => {
  it('returns value for a simple key', () => {
    expect(getNestedValue({name: 'Chem'}, 'name')).toBe('Chem');
  });

  it('returns value for a nested path', () => {
    expect(getNestedValue({a: {b: {c: 42}}}, 'a.b.c')).toBe(42);
  });

  it('returns undefined for a missing key', () => {
    expect(getNestedValue({name: 'Chem'}, 'version')).toBeUndefined();
  });

  it('returns undefined when a mid-path segment is null', () => {
    expect(getNestedValue({a: null}, 'a.b')).toBeUndefined();
  });

  it('returns undefined when a mid-path segment is missing', () => {
    expect(getNestedValue({a: {}}, 'a.b.c')).toBeUndefined();
  });

  it('returns undefined for an empty path (splits to empty string key)', () => {
    expect(getNestedValue({x: 1}, '')).toBeUndefined();
  });
});

const pkg = (overrides: Record<string, any>) => ({
  dir: '/tmp/pkg',
  name: overrides.name ?? 'test-pkg',
  friendlyName: overrides.friendlyName ?? overrides.name ?? 'Test Pkg',
  version: overrides.version ?? '1.0.0',
  packageJson: overrides,
});

describe('applyFilter', () => {
  const packages = [
    pkg({name: 'Chem', version: '1.5.0', category: 'Cheminformatics'}),
    pkg({name: 'Bio', version: '2.0.0', category: 'Bioinformatics'}),
    pkg({name: 'PowerGrid', version: '1.5.0', category: 'Viewers'}),
  ];

  it('returns all packages when filter matches all', () => {
    expect(applyFilter(packages, 'name:.')).toHaveLength(3);
  });

  it('filters by exact name match', () => {
    const result = applyFilter(packages, 'name:^Chem$');
    expect(result).toHaveLength(1);
    expect(result[0].name).toBe('Chem');
  });

  it('filters by partial name (regex substring)', () => {
    const result = applyFilter(packages, 'name:Bio');
    expect(result).toHaveLength(1);
    expect(result[0].name).toBe('Bio');
  });

  it('returns empty array when nothing matches', () => {
    expect(applyFilter(packages, 'name:NOMATCH')).toHaveLength(0);
  });

  it('filters by version', () => {
    const result = applyFilter(packages, 'version:^1\\.5');
    expect(result).toHaveLength(2);
    expect(result.map((p) => p.name)).toEqual(expect.arrayContaining(['Chem', 'PowerGrid']));
  });

  it('applies && conjunction (both conditions must match)', () => {
    const result = applyFilter(packages, 'name:Chem && version:1\\.5');
    expect(result).toHaveLength(1);
    expect(result[0].name).toBe('Chem');
  });

  it('returns empty when one part of && conjunction fails', () => {
    expect(applyFilter(packages, 'name:Chem && version:^2')).toHaveLength(0);
  });

  it('filters by nested field', () => {
    const withNested = [
      pkg({name: 'A', datagrok: {apiVersion: '1.0'}}),
      pkg({name: 'B', datagrok: {apiVersion: '2.0'}}),
    ];
    const result = applyFilter(withNested, 'datagrok.apiVersion:^1');
    expect(result).toHaveLength(1);
    expect(result[0].name).toBe('A');
  });

  it('returns empty when field does not exist', () => {
    expect(applyFilter(packages, 'nonexistent:anything')).toHaveLength(0);
  });

  it('treats filter with no colon as field name with match-all pattern', () => {
    // No colon → field = whole string, pattern = /./  (matches any value)
    // The function returns packages where the field exists and is non-empty
    const result = applyFilter(packages, 'name');
    expect(result).toHaveLength(3);
  });
});
