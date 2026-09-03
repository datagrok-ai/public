/* Binding discovery: every module under a binding source is imported for its side effects
   (registrations) and inspected for exported step functions, so the compiler can emit a named
   import for each step it uses and a side-effect import for each module that registers nouns.
   Library bindings are the compiled `.js` of this build (or `.ts` from source); a project's own
   bindings are TypeScript, loaded through tsx. */
import {existsSync, readdirSync, statSync} from 'node:fs';
import {join} from 'node:path';
import {pathToFileURL} from 'node:url';
import {BindingSource, LIB_DIR, LIB_ROOT} from './project.js';
import {allDatasets, allElements, allKinds, parameterTypes, StepDef, StepFn, steps} from './registry.js';

export interface BindingModule {
  file: string;
  /** What a generated spec imports: a package subpath, or the file itself for project bindings. */
  specifier: string;
  exports: Record<string, unknown>;
  stepDefs: StepDef[];
  /** Registered elements, contexts, kinds, datasets or parameter types at import. */
  registers: boolean;
}

export interface Bindings {
  modules: BindingModule[];
  exportOf: Map<StepFn, {module: BindingModule; name: string}>;
}

const EXTENSIONS = ['.ts', '.mts', '.js', '.mjs'];

export function listFiles(dir: string, ext: string | string[] = EXTENSIONS): string[] {
  if (!existsSync(dir))
    return [];
  const exts = Array.isArray(ext) ? ext : [ext];
  const out: string[] = [];
  for (const name of readdirSync(dir).sort()) {
    const path = join(dir, name);
    if (statSync(path).isDirectory()) {
      if (name !== 'node_modules' && !name.startsWith('.'))
        out.push(...listFiles(path, exts));
    }
    else if (exts.some((e) => name.endsWith(e)) && !/\.d\.[cm]?ts$/.test(name))
      out.push(path);
  }
  return out;
}

let tsReady = LIB_DIR === LIB_ROOT;

async function ensureTypeScript(): Promise<void> {
  if (tsReady)
    return;
  const {register} = await import('tsx/esm/api');
  register();
  tsReady = true;
}

export async function loadBindings(sources: BindingSource[]): Promise<Bindings> {
  const modules: BindingModule[] = [];
  const exportOf = new Map<StepFn, {module: BindingModule; name: string}>();
  for (const source of sources) {
    for (const file of listFiles(source.dir)) {
      if (/\.[cm]?ts$/.test(file))
        await ensureTypeScript();
      const before = {steps: steps.length, elements: allElements().length, kinds: allKinds().length,
        datasets: allDatasets().length, params: parameterTypes.length};
      const exports = await import(pathToFileURL(file).href) as Record<string, unknown>;
      const stepDefs = steps.slice(before.steps);
      const registers = allElements().length !== before.elements || allKinds().length !== before.kinds ||
        allDatasets().length !== before.datasets || parameterTypes.length !== before.params;
      const module: BindingModule = {file, specifier: source.specifierFor(file), exports, stepDefs, registers};
      modules.push(module);
      for (const [name, value] of Object.entries(exports)) {
        for (const def of stepDefs) {
          if (def.fn === value)
            exportOf.set(def.fn, {module, name});
        }
      }
    }
  }
  return {modules, exportOf};
}
