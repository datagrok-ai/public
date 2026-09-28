import {describe, it, expect, beforeEach, afterEach, vi} from 'vitest';
import * as os from 'os';
import * as fs from 'fs';
import * as path from 'path';
import {parseFuncCall, buildInlineManifest, resolveManifestSources} from '../commands/server';

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

describe('buildInlineManifest', () => {
  it('sets action to entity.verb', () => {
    const result = buildInlineManifest('users', 'delete', ['abc']);
    expect(result.operations[0].action).toBe('users.delete');
  });

  it('maps string args to {id} param for non-file entities', () => {
    const result = buildInlineManifest('users', 'delete', ['id1', 'id2']);
    expect(result.operations).toEqual([
      {id: 'op0', action: 'users.delete', params: {id: 'id1'}},
      {id: 'op1', action: 'users.delete', params: {id: 'id2'}},
    ]);
  });

  it('maps string args to {path} param for the files entity', () => {
    const result = buildInlineManifest('files', 'delete', ['System:AppData/a.txt', 'System:AppData/b.txt']);
    expect(result.operations).toEqual([
      {id: 'op0', action: 'files.delete', params: {path: 'System:AppData/a.txt'}},
      {id: 'op1', action: 'files.delete', params: {path: 'System:AppData/b.txt'}},
    ]);
  });

  it('passes object array args through as params directly', () => {
    const objs = [{name: 'Alice', email: 'a@x.com'}, {name: 'Bob', email: 'b@x.com'}];
    const result = buildInlineManifest('users', 'create', objs);
    expect(result.operations).toEqual([
      {id: 'op0', action: 'users.create', params: {name: 'Alice', email: 'a@x.com'}},
      {id: 'op1', action: 'users.create', params: {name: 'Bob', email: 'b@x.com'}},
    ]);
  });

  it('assigns sequential op ids starting at op0', () => {
    const result = buildInlineManifest('groups', 'delete', ['x', 'y', 'z']);
    expect(result.operations.map((o) => o.id)).toEqual(['op0', 'op1', 'op2']);
  });

  it('returns a single operation for a single arg', () => {
    const result = buildInlineManifest('connections', 'delete', ['conn-id']);
    expect(result.operations).toHaveLength(1);
  });
});

describe('resolveManifestSources', () => {
  let tmpFile: string;

  beforeEach(() => {
    tmpFile = path.join(os.tmpdir(), `grok-batch-test-${Date.now()}.bin`);
  });

  afterEach(() => {
    if (fs.existsSync(tmpFile)) fs.unlinkSync(tmpFile);
  });

  it('passes through a manifest with no files.put operations unchanged', () => {
    const manifest = {
      operations: [{id: 'op0', action: 'users.delete', params: {id: 'abc'}}],
    };
    expect(resolveManifestSources(manifest)).toEqual(manifest);
  });

  it('passes through a files.put operation that has no source field', () => {
    const manifest = {
      operations: [{id: 'op0', action: 'files.put', params: {path: 'System:AppData/f.txt', content: 'aGk='}}],
    };
    expect(resolveManifestSources(manifest)).toEqual(manifest);
  });

  it('reads the source file, base64-encodes its contents, and removes the source key', () => {
    const data = 'hello world';
    fs.writeFileSync(tmpFile, data);

    const manifest = {
      operations: [{id: 'op0', action: 'files.put', params: {path: 'System:AppData/f.txt', source: tmpFile}}],
    };
    const result = resolveManifestSources(manifest);
    const params = result.operations[0].params as any;
    expect(params.content).toBe(Buffer.from(data).toString('base64'));
    expect(params.source).toBeUndefined();
  });

  it('preserves other params alongside the injected content', () => {
    fs.writeFileSync(tmpFile, 'data');

    const manifest = {
      operations: [{id: 'op0', action: 'files.put', params: {path: 'System:AppData/f.txt', source: tmpFile, extra: 'val'}}],
    };
    const result = resolveManifestSources(manifest);
    const params = result.operations[0].params as any;
    expect(params.path).toBe('System:AppData/f.txt');
    expect(params.extra).toBe('val');
  });

  it('leaves non-files.put operations untouched in a mixed manifest', () => {
    fs.writeFileSync(tmpFile, 'x');

    const manifest = {
      operations: [
        {id: 'op0', action: 'users.delete', params: {id: 'u1'}},
        {id: 'op1', action: 'files.put', params: {path: 'System:AppData/f.txt', source: tmpFile}},
      ],
    };
    const result = resolveManifestSources(manifest);
    expect(result.operations[0]).toEqual({id: 'op0', action: 'users.delete', params: {id: 'u1'}});
    expect((result.operations[1].params as any).source).toBeUndefined();
    expect((result.operations[1].params as any).content).toBeDefined();
  });

  it('correctly encodes binary file contents', () => {
    const bytes = Buffer.from([0x00, 0xff, 0x10, 0xab]);
    fs.writeFileSync(tmpFile, bytes);

    const manifest = {
      operations: [{id: 'op0', action: 'files.put', params: {path: 'System:AppData/f.bin', source: tmpFile}}],
    };
    const result = resolveManifestSources(manifest);
    const decoded = Buffer.from((result.operations[0].params as any).content, 'base64');
    expect(decoded).toEqual(bytes);
  });
});

describe('observability dispatch', () => {
  it('answers `grok s alerts --help` with the command usage', async () => {
    const {server} = await import('../commands/server');
    const {ALERTS_USAGE} = await import('../commands/server-alerts');
    const lines: string[] = [];
    const log = vi.spyOn(console, 'log').mockImplementation((s: any) => { lines.push(String(s)); });
    try {
      expect(await server({_: ['s', 'alerts'], help: true})).toBe(true);
    } finally {
      log.mockRestore();
    }
    expect(lines).toEqual([ALERTS_USAGE]);
  });

  it('answers a bare `grok s --help` with success', async () => {
    const {server, HELP_SERVER} = await import('../commands/server');
    const log = vi.spyOn(console, 'log').mockImplementation(() => {});
    try {
      expect(await server({_: ['s'], help: true})).toBe(true);
      expect(log).toHaveBeenCalledWith(HELP_SERVER);
    } finally {
      log.mockRestore();
    }
  });

  it('refuses a repeated --host outside alerts, errors and logger', async () => {
    const {server} = await import('../commands/server');
    const err: string[] = [];
    const write = vi.spyOn(process.stderr, 'write').mockImplementation((s: any) => { err.push(String(s)); return true; });
    const before = process.exitCode;
    try {
      expect(await server({_: ['s', 'users', 'list'], host: ['a', 'b']})).toBe(true);
      expect(process.exitCode).toBe(1);
    } finally {
      write.mockRestore();
      process.exitCode = before;
    }
    expect(err.join('')).toMatch(/--host may repeat only for alerts, errors and logger/);
  });
});
