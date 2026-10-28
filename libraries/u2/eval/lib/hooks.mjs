// Loader hooks: u2 from its test build (.test-build), .ts transpiled as u2 compiles, .css installed as package CSS.
import {existsSync, readFileSync} from 'node:fs';
import {fileURLToPath, pathToFileURL} from 'node:url';
import ts from 'typescript';

const U2 = new URL('../../', import.meta.url);
const SRC = fileURLToPath(new URL('src/', U2));
const BUILT = fileURLToPath(new URL('.test-build/src/', U2));
const PLATFORM = new URL('./platform.mjs', import.meta.url).href;
/** The page's platform modules: the stub's, plus what AGENT-PROMPT.md promises and the stub lacks. */
const WRAPPERS = {
  'datagrok-api/dg': ['eval:dg', 'export * from \'datagrok-api/dg\';\n' +
    `export {TableView, View, ViewBase} from '${PLATFORM}';\n`],
  'datagrok-api/grok': ['eval:grok', 'export * from \'datagrok-api/grok\';\n' +
    'export const chem = globalThis.__evalChem ??= {};\n'],
  'datagrok-api/ui': ['eval:ui', 'export * from \'datagrok-api/ui\';\n' +
    'export const pickFileFromFiles = () => globalThis.grok_UI_PickFileFromFiles();\n'],
};
const SOURCES = new Map(Object.values(WRAPPERS));

const COMPILER_OPTIONS = {
  target: ts.ScriptTarget.ES2022, module: ts.ModuleKind.ESNext, useDefineForClassFields: false,
  experimentalDecorators: true, esModuleInterop: true,
};

/** u2's entry points (`exports` in its package.json) as the source folders they stand for. */
const ENTRIES = {'': 'src', '/dg': 'src/dg', '/filter': 'src/core/filter', '/dg/filter': 'src/dg/filter',
  '/testing': 'src/testing'};

function u2Url(subpath) {
  let rest = ENTRIES[subpath] ?? subpath.replace(/^\//, '').replace(/\.ts$/, '.js');
  if (rest.startsWith('src'))
    rest = '.test-build/' + rest;
  if (!rest.endsWith('.js'))
    rest += existsSync(fileURLToPath(new URL(rest, U2))) ? '/index.js' : '.js';
  return new URL(rest, U2).href;
}

export async function resolve(specifier, context, next) {
  if (specifier === '@datagrok-libraries/u2' || specifier.startsWith('@datagrok-libraries/u2/'))
    return {url: u2Url(specifier.slice('@datagrok-libraries/u2'.length)), format: 'module', shortCircuit: true};
  const wrapper = WRAPPERS[specifier]?.[0];
  if (wrapper !== undefined && context.parentURL !== wrapper)
    return {url: wrapper, format: 'module', shortCircuit: true};
  if (/^\.\.?\//.test(specifier) && context.parentURL?.startsWith('file:')) {
    const url = new URL(specifier, context.parentURL);
    if (specifier.endsWith('.css'))
      return {url: url.href, format: 'module', shortCircuit: true};
    const path = fileURLToPath(url);
    if (path.startsWith(SRC))
      return {url: pathToFileURL(BUILT + path.slice(SRC.length).replace(/\.ts$/, '.js')).href, format: 'module',
        shortCircuit: true};
    for (const candidate of [path, path.replace(/\.js$/, '.ts'), path + '.ts']) {
      if (candidate.endsWith('.ts') && existsSync(candidate))
        return {url: pathToFileURL(candidate).href, format: 'module', shortCircuit: true};
    }
  }
  return next(specifier, context);
}

export async function load(url, context, next) {
  if (SOURCES.has(url))
    return {format: 'module', shortCircuit: true, source: SOURCES.get(url)};
  if (url.startsWith('file:') && url.endsWith('.ts')) {
    const text = readFileSync(fileURLToPath(url), 'utf8');
    const out = ts.transpileModule(text, {compilerOptions: COMPILER_OPTIONS, fileName: fileURLToPath(url)});
    return {format: 'module', shortCircuit: true, source: out.outputText};
  }
  if (url.startsWith('file:') && url.endsWith('.css')) {
    const css = JSON.stringify(readFileSync(fileURLToPath(url), 'utf8'));
    return {format: 'module', shortCircuit: true, source:
      `const s = document.createElement('style');\ns.setAttribute('data-package-css', ${JSON.stringify(url)});\n` +
      `s.textContent = ${css};\ndocument.head.append(s);\nexport default {};\n`};
  }
  return next(url, context);
}
