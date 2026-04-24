"use strict";

var _vitest = require("vitest");
var os = _interopRequireWildcard(require("os"));
var fs = _interopRequireWildcard(require("fs"));
var path = _interopRequireWildcard(require("path"));
var _server = require("../commands/server");
function _interopRequireWildcard(e, t) { if ("function" == typeof WeakMap) var r = new WeakMap(), n = new WeakMap(); return (_interopRequireWildcard = function (e, t) { if (!t && e && e.__esModule) return e; var o, i, f = { __proto__: null, default: e }; if (null === e || "object" != typeof e && "function" != typeof e) return f; if (o = t ? n : r) { if (o.has(e)) return o.get(e); o.set(e, f); } for (const t in e) "default" !== t && {}.hasOwnProperty.call(e, t) && ((i = (o = Object.defineProperty) && Object.getOwnPropertyDescriptor(e, t)) && (i.get || i.set) ? o(f, t, i) : f[t] = e[t]); return f; })(e, t); }
(0, _vitest.describe)('parseFuncCall', () => {
  (0, _vitest.it)('parses a single string argument', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Chem:smilesToMw("ccc")')).toEqual({
      name: 'Chem:smilesToMw',
      params: {
        '0': 'ccc'
      }
    });
  });
  (0, _vitest.it)('parses a single-quoted string argument', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)("Pkg:fn('hello')")).toEqual({
      name: 'Pkg:fn',
      params: {
        '0': 'hello'
      }
    });
  });
  (0, _vitest.it)('parses a numeric argument', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn(42)')).toEqual({
      name: 'Pkg:fn',
      params: {
        '0': 42
      }
    });
  });
  (0, _vitest.it)('parses multiple positional arguments', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn(1, "hello", 3.14)')).toEqual({
      name: 'Pkg:fn',
      params: {
        '0': 1,
        '1': 'hello',
        '2': 3.14
      }
    });
  });
  (0, _vitest.it)('parses an object argument with quoted keys', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn({"a":5,"b":22})')).toEqual({
      name: 'Pkg:fn',
      params: {
        a: 5,
        b: 22
      }
    });
  });
  (0, _vitest.it)('parses an object argument with unquoted keys', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn({a:5,b:22})')).toEqual({
      name: 'Pkg:fn',
      params: {
        a: 5,
        b: 22
      }
    });
  });
  (0, _vitest.it)('parses an object argument with a boolean value', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn({unquoted:true})')).toEqual({
      name: 'Pkg:fn',
      params: {
        unquoted: true
      }
    });
  });
  (0, _vitest.it)('returns empty params for a no-argument call', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn()')).toEqual({
      name: 'Pkg:fn',
      params: {}
    });
  });
  (0, _vitest.it)('returns empty params when there are no parentheses', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('Pkg:fn')).toEqual({
      name: 'Pkg:fn',
      params: {}
    });
  });
  (0, _vitest.it)('handles a package-less function name', () => {
    (0, _vitest.expect)((0, _server.parseFuncCall)('myFunc("arg")')).toEqual({
      name: 'myFunc',
      params: {
        '0': 'arg'
      }
    });
  });
});
(0, _vitest.describe)('buildInlineManifest', () => {
  (0, _vitest.it)('sets action to entity.verb', () => {
    const result = (0, _server.buildInlineManifest)('users', 'delete', ['abc']);
    (0, _vitest.expect)(result.operations[0].action).toBe('users.delete');
  });
  (0, _vitest.it)('maps string args to {id} param for non-file entities', () => {
    const result = (0, _server.buildInlineManifest)('users', 'delete', ['id1', 'id2']);
    (0, _vitest.expect)(result.operations).toEqual([{
      id: 'op0',
      action: 'users.delete',
      params: {
        id: 'id1'
      }
    }, {
      id: 'op1',
      action: 'users.delete',
      params: {
        id: 'id2'
      }
    }]);
  });
  (0, _vitest.it)('maps string args to {path} param for the files entity', () => {
    const result = (0, _server.buildInlineManifest)('files', 'delete', ['System:AppData/a.txt', 'System:AppData/b.txt']);
    (0, _vitest.expect)(result.operations).toEqual([{
      id: 'op0',
      action: 'files.delete',
      params: {
        path: 'System:AppData/a.txt'
      }
    }, {
      id: 'op1',
      action: 'files.delete',
      params: {
        path: 'System:AppData/b.txt'
      }
    }]);
  });
  (0, _vitest.it)('passes object array args through as params directly', () => {
    const objs = [{
      name: 'Alice',
      email: 'a@x.com'
    }, {
      name: 'Bob',
      email: 'b@x.com'
    }];
    const result = (0, _server.buildInlineManifest)('users', 'create', objs);
    (0, _vitest.expect)(result.operations).toEqual([{
      id: 'op0',
      action: 'users.create',
      params: {
        name: 'Alice',
        email: 'a@x.com'
      }
    }, {
      id: 'op1',
      action: 'users.create',
      params: {
        name: 'Bob',
        email: 'b@x.com'
      }
    }]);
  });
  (0, _vitest.it)('assigns sequential op ids starting at op0', () => {
    const result = (0, _server.buildInlineManifest)('groups', 'delete', ['x', 'y', 'z']);
    (0, _vitest.expect)(result.operations.map(o => o.id)).toEqual(['op0', 'op1', 'op2']);
  });
  (0, _vitest.it)('returns a single operation for a single arg', () => {
    const result = (0, _server.buildInlineManifest)('connections', 'delete', ['conn-id']);
    (0, _vitest.expect)(result.operations).toHaveLength(1);
  });
});
(0, _vitest.describe)('resolveManifestSources', () => {
  let tmpFile;
  (0, _vitest.beforeEach)(() => {
    tmpFile = path.join(os.tmpdir(), `grok-batch-test-${Date.now()}.bin`);
  });
  (0, _vitest.afterEach)(() => {
    if (fs.existsSync(tmpFile)) fs.unlinkSync(tmpFile);
  });
  (0, _vitest.it)('passes through a manifest with no files.put operations unchanged', () => {
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'users.delete',
        params: {
          id: 'abc'
        }
      }]
    };
    (0, _vitest.expect)((0, _server.resolveManifestSources)(manifest)).toEqual(manifest);
  });
  (0, _vitest.it)('passes through a files.put operation that has no source field', () => {
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'files.put',
        params: {
          path: 'System:AppData/f.txt',
          content: 'aGk='
        }
      }]
    };
    (0, _vitest.expect)((0, _server.resolveManifestSources)(manifest)).toEqual(manifest);
  });
  (0, _vitest.it)('reads the source file, base64-encodes its contents, and removes the source key', () => {
    const data = 'hello world';
    fs.writeFileSync(tmpFile, data);
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'files.put',
        params: {
          path: 'System:AppData/f.txt',
          source: tmpFile
        }
      }]
    };
    const result = (0, _server.resolveManifestSources)(manifest);
    const params = result.operations[0].params;
    (0, _vitest.expect)(params.content).toBe(Buffer.from(data).toString('base64'));
    (0, _vitest.expect)(params.source).toBeUndefined();
  });
  (0, _vitest.it)('preserves other params alongside the injected content', () => {
    fs.writeFileSync(tmpFile, 'data');
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'files.put',
        params: {
          path: 'System:AppData/f.txt',
          source: tmpFile,
          extra: 'val'
        }
      }]
    };
    const result = (0, _server.resolveManifestSources)(manifest);
    const params = result.operations[0].params;
    (0, _vitest.expect)(params.path).toBe('System:AppData/f.txt');
    (0, _vitest.expect)(params.extra).toBe('val');
  });
  (0, _vitest.it)('leaves non-files.put operations untouched in a mixed manifest', () => {
    fs.writeFileSync(tmpFile, 'x');
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'users.delete',
        params: {
          id: 'u1'
        }
      }, {
        id: 'op1',
        action: 'files.put',
        params: {
          path: 'System:AppData/f.txt',
          source: tmpFile
        }
      }]
    };
    const result = (0, _server.resolveManifestSources)(manifest);
    (0, _vitest.expect)(result.operations[0]).toEqual({
      id: 'op0',
      action: 'users.delete',
      params: {
        id: 'u1'
      }
    });
    (0, _vitest.expect)(result.operations[1].params.source).toBeUndefined();
    (0, _vitest.expect)(result.operations[1].params.content).toBeDefined();
  });
  (0, _vitest.it)('correctly encodes binary file contents', () => {
    const bytes = Buffer.from([0x00, 0xff, 0x10, 0xab]);
    fs.writeFileSync(tmpFile, bytes);
    const manifest = {
      operations: [{
        id: 'op0',
        action: 'files.put',
        params: {
          path: 'System:AppData/f.bin',
          source: tmpFile
        }
      }]
    };
    const result = (0, _server.resolveManifestSources)(manifest);
    const decoded = Buffer.from(result.operations[0].params.content, 'base64');
    (0, _vitest.expect)(decoded).toEqual(bytes);
  });
});