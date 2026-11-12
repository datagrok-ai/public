// One TypeScript program over all solutions, typed against u2's sources; u2's own files are not checked.
import {fileURLToPath} from 'node:url';
import ts from 'typescript';

const U2 = fileURLToPath(new URL('../../', import.meta.url));

const OPTIONS = {
  target: ts.ScriptTarget.ES2022, module: ts.ModuleKind.ESNext, moduleResolution: ts.ModuleResolutionKind.Bundler,
  lib: ['lib.esnext.d.ts', 'lib.dom.d.ts', 'lib.dom.iterable.d.ts'], strict: true, useDefineForClassFields: false,
  experimentalDecorators: true, esModuleInterop: true, skipLibCheck: true, noEmit: true, types: [],
  baseUrl: U2, paths: {
    '@datagrok-libraries/u2': ['src/index.ts'],
    '@datagrok-libraries/u2/dg': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/filter': ['src/core/filter/index.ts'],
    '@datagrok-libraries/u2/dg/filter': ['src/dg/filter/index.ts'],
    '@datagrok-libraries/u2/testing': ['src/testing/index.ts'],
    '@datagrok-libraries/u2/src/dg/index.js': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/src/dg': ['src/dg/index.ts'],
    '@datagrok-libraries/u2/*': ['*'],
    'datagrok-api/dg': ['node_modules/datagrok-api/dist/dg.d.ts'],
    'datagrok-api/grok': ['node_modules/datagrok-api/dist/grok.d.ts'],
    'datagrok-api/ui': ['node_modules/datagrok-api/dist/ui.d.ts'],
    'datagrok-api/u2core': ['node_modules/datagrok-api/dist/src/u2core/index.d.ts'],
  },
};

const CSS_MODULE = 'declare module \'*.css\' { const css: unknown; export default css; }\n';

export function typecheck(files) {
  const host = ts.createCompilerHost(OPTIONS);
  const shim = U2 + 'eval/lib/__css.d.ts';
  const read = host.readFile.bind(host);
  host.readFile = (f) => f === shim ? CSS_MODULE : read(f);
  const exists = host.fileExists.bind(host);
  host.fileExists = (f) => f === shim || exists(f);
  const getSource = host.getSourceFile.bind(host);
  host.getSourceFile = (f, v, onError, create) =>
    f === shim ? ts.createSourceFile(f, CSS_MODULE, v) : getSource(f, v, onError, create);
  const program = ts.createProgram([...files, shim], OPTIONS, host);
  const result = new Map();
  for (const file of files) {
    const sf = program.getSourceFile(file);
    const diagnostics = sf === undefined ? [] :
      [...program.getSyntacticDiagnostics(sf), ...program.getSemanticDiagnostics(sf)];
    result.set(file, diagnostics.map((d) => {
      const line = d.start === undefined ? '?' : sf.getLineAndCharacterOfPosition(d.start).line + 1;
      return `${line}: ${ts.flattenDiagnosticMessageText(d.messageText, ' ')}`;
    }));
  }
  return result;
}
