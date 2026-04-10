import {describe, it, expect} from 'vitest';
import {parseFuncCall} from '../commands/server';

describe('parseFuncCall', () => {
  it('parses a single string argument', () => {
    expect(parseFuncCall('Chem:smilesToMw("ccc")')).toEqual({
      name: 'Chem:smilesToMw',
      params: {'0': 'ccc'},
    });
  });

  it('parses a single-quoted string argument', () => {
    expect(parseFuncCall("Pkg:fn('hello')")).toEqual({
      name: 'Pkg:fn',
      params: {'0': 'hello'},
    });
  });

  it('parses a numeric argument', () => {
    expect(parseFuncCall('Pkg:fn(42)')).toEqual({
      name: 'Pkg:fn',
      params: {'0': 42},
    });
  });

  it('parses multiple positional arguments', () => {
    expect(parseFuncCall('Pkg:fn(1, "hello", 3.14)')).toEqual({
      name: 'Pkg:fn',
      params: {'0': 1, '1': 'hello', '2': 3.14},
    });
  });

  it('parses an object argument with quoted keys', () => {
    expect(parseFuncCall('Pkg:fn({"a":5,"b":22})')).toEqual({
      name: 'Pkg:fn',
      params: {a: 5, b: 22},
    });
  });

  it('parses an object argument with unquoted keys', () => {
    expect(parseFuncCall('Pkg:fn({a:5,b:22})')).toEqual({
      name: 'Pkg:fn',
      params: {a: 5, b: 22},
    });
  });

  it('parses an object argument with a boolean value', () => {
    expect(parseFuncCall('Pkg:fn({unquoted:true})')).toEqual({
      name: 'Pkg:fn',
      params: {unquoted: true},
    });
  });

  it('returns empty params for a no-argument call', () => {
    expect(parseFuncCall('Pkg:fn()')).toEqual({
      name: 'Pkg:fn',
      params: {},
    });
  });

  it('returns empty params when there are no parentheses', () => {
    expect(parseFuncCall('Pkg:fn')).toEqual({
      name: 'Pkg:fn',
      params: {},
    });
  });

  it('handles a package-less function name', () => {
    expect(parseFuncCall('myFunc("arg")')).toEqual({
      name: 'myFunc',
      params: {'0': 'arg'},
    });
  });
});
