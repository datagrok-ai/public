/* Feature model → one Playwright spec. Deterministic: the output is a pure function of the feature
   text and the loaded bindings; element phrases (`{element}`, `{widget}`) are emitted as `el('…')` and datasets as `ds('…')`
   (names, never selectors), so a registry fix never forces a regeneration. A feature's scenarios
   share one browser page through `feature(test)` (see runtime/harness.ts): Playwright still runs
   one test per scenario. A step declared with `enters` switches the vocabulary: the compiler
   validates the following phrases against that context and emits `enter(page, '…')` so the
   runtime resolves them the same way. */
import {dirname, isAbsolute, join, relative, sep} from 'node:path';
import type {Bindings} from './discover.js';
import type {FeatureModel, ScenarioModel, StepModel} from './gherkin.js';
import type {MatchedArg, StepMatcher} from './match.js';
import {describeNoun, NounError, parseNoun, usesKinds} from './nouns.js';
import {PACKAGE_NAME} from './project.js';
import {ContextEntry, lookupContext, lookupDataset} from './registry.js';

export type DiagnosticLevel = 'error' | 'warning' | 'note';

export interface Diagnostic {
  file: string;
  line: number;
  level: DiagnosticLevel;
  message: string;
}

export interface CompileContext {
  /** The project root: `features/` and `generated/` are right under it. */
  root: string;
  matcher: StepMatcher;
  bindings: Bindings;
}

export interface CompiledFeature {
  feature: FeatureModel;
  outFile: string;
  code: string;
  diagnostics: Diagnostic[];
}

/** A feature whose scenarios run in order on one shell state, the Background once (see `journey`). */
export const JOURNEY_TAG = '@journey';
export const GENERATED_DIR = 'generated';
export const FEATURES_DIR = 'features';
export const RUNTIME_SPECIFIER = `${PACKAGE_NAME}/runtime`;

function posix(p: string): string {
  return p.split(sep).join('/');
}

/** `<root>/features/x/y.feature` → `<root>/generated/x/y.test.ts`. */
export function outFileFor(root: string, featurePath: string): string {
  const rel = relative(join(root, FEATURES_DIR), featurePath);
  return join(root, GENERATED_DIR, rel.replace(/\.feature$/, '.test.ts'));
}

/** A package subpath as is; a project file relative to the spec, with the extension a module
 * resolver expects (`.js` for `.ts` — Playwright maps it back). */
function importSpecifier(from: string, specifier: string): string {
  if (!isAbsolute(specifier))
    return specifier;
  const p = posix(relative(dirname(from), specifier)).replace(/\.mts$/, '.mjs').replace(/\.ts$/, '.js');
  return p.startsWith('.') ? p : './' + p;
}

const IDENT = /^[A-Za-z_$][\w$]*$/;

interface ScenarioState {
  context?: ContextEntry;
}

export function compileFeature(feature: FeatureModel, ctx: CompileContext): CompiledFeature {
  const outFile = outFileFor(ctx.root, feature.path);
  const diagnostics: Diagnostic[] = [];
  const named = new Map<string, Set<string>>();
  const helpers = new Set<string>(['feature']);
  const relPath = posix(relative(ctx.root, feature.path));
  const diag = (line: number, level: DiagnosticLevel, message: string) =>
    diagnostics.push({file: relPath, line, level, message});

  const emitArg = (arg: MatchedArg, step: StepModel, context: ContextEntry | undefined): string => {
    switch (arg.type) {
      case 'element':
      case 'widget': {
        helpers.add('el');
        const phrase = String(arg.value);
        try {
          const ref = parseNoun(phrase, context);
          if (usesKinds(ref))
            diag(step.line, 'note', `"${phrase}" → ${describeNoun(ref)}`);
        } catch (e) {
          if (e instanceof NounError)
            diag(step.line, 'error', `element "${phrase}": ${e.message}`);
          else
            throw e;
        }
        return `el(${JSON.stringify(phrase)})`;
      }
      case 'dataset': {
        helpers.add('ds');
        const name = String(arg.value);
        if (!lookupDataset(name))
          diag(step.line, 'error', `dataset "${name}" is not registered`);
        return `ds(${JSON.stringify(name)})`;
      }
      case 'int':
      case 'float':
      case 'double':
        return String(arg.value);
      default:
        return JSON.stringify(arg.value);
    }
  };

  const emitStep = (step: StepModel, indent: string, state: ScenarioState): string[] => {
    const title = JSON.stringify(`${step.keyword} ${step.text}`);
    const result = ctx.matcher.match(step.text);
    if (result.ambiguous) {
      diag(step.line, 'error', `ambiguous step "${step.text}": ` +
        result.ambiguous.map((m) => `"${m.def.expression}"`).join(' | '));
      return [`${indent}await test.step(${title}, () => { throw new Error('ambiguous step'); });`];
    }
    if (!result.match) {
      const near = ctx.matcher.suggest(step.text).map((d) => `"${d.expression}"`).join(', ');
      diag(step.line, 'error', `no step definition matches "${step.text}"` + (near ? ` — nearest: ${near}` : ''));
      return [`${indent}await test.step(${title}, () => { throw new Error('no step definition matches this step'); });`];
    }
    const {def, args} = result.match;
    const exported = ctx.bindings.exportOf.get(def.fn);
    if (!exported || !IDENT.test(exported.name)) {
      diag(step.line, 'error', `the definition of "${def.expression}" must be an exported const of its module`);
      return [`${indent}await test.step(${title}, () => { throw new Error('step definition is not exported'); });`];
    }
    if (!named.has(exported.module.specifier))
      named.set(exported.module.specifier, new Set());
    named.get(exported.module.specifier)!.add(exported.name);
    const call = [...args.map((a) => emitArg(a, step, state.context))];
    if (step.table)
      call.push(JSON.stringify(step.table));
    if (step.docString !== undefined)
      call.push(JSON.stringify(step.docString));
    const out = [`${indent}await test.step(${title}, () => ${exported.name}(page${call.map((c) => ', ' + c).join('')}));`];
    if (def.meta.enters !== undefined) {
      const entered = lookupContext(def.meta.enters);
      if (!entered)
        diag(step.line, 'error', `"${def.expression}" enters "${def.meta.enters}", which is not a registered context`);
      else {
        state.context = entered;
        helpers.add('enter');
        out.push(`${indent}enter(page, ${JSON.stringify(def.meta.enters)});`);
      }
    }
    else if (def.meta.leaves && state.context) {
      state.context = undefined;
      helpers.add('leave');
      out.push(`${indent}leave(page);`);
    }
    return out;
  };

  const seen = new Map<string, number>();
  const uniqueTitle = (name: string): string => {
    const n = (seen.get(name) ?? 0) + 1;
    seen.set(name, n);
    return n === 1 ? name : `${name} (${n})`;
  };

  const body: string[] = feature.tags.includes(JOURNEY_TAG) ?
    emitJourney(feature, uniqueTitle, emitStep, helpers) :
    feature.scenarios.flatMap((scenario) => emitScenario(scenario, feature, uniqueTitle, emitStep));

  const lines: string[] = [];
  const realizes = [...feature.tags, ...feature.scenarios.flatMap((s) => s.tags)]
    .filter((t) => t.startsWith('@realizes:')).map((t) => t.slice('@realizes:'.length));
  lines.push('/* eslint-disable max-len */', '/* eslint-disable comma-spacing */', '/* eslint-disable quotes */');
  lines.push('/* ---');
  lines.push(`generated: ${relPath}`);
  lines.push(`generator: ${PACKAGE_NAME} — do not edit; run \`grok-bdd compile\` to regenerate`);
  if (realizes.length > 0)
    lines.push(`sub_features_covered: [${[...new Set(realizes)].join(', ')}]`);
  lines.push('--- */');
  lines.push(`import {test} from '@playwright/test';`);
  const sideEffects = ctx.bindings.modules.filter((m) => m.registers && !named.has(m.specifier))
    .map((m) => importSpecifier(outFile, m.specifier)).sort();
  for (const specifier of sideEffects)
    lines.push(`import '${specifier}';`);
  for (const [specifier, names] of [...named.entries()].map(([s, n]) => [importSpecifier(outFile, s), n] as const).sort((a, b) => a[0].localeCompare(b[0])))
    lines.push(`import {${[...names].sort().join(', ')}} from '${specifier}';`);
  lines.push(`import {${[...helpers].sort().join(', ')}} from '${RUNTIME_SPECIFIER}';`);
  lines.push('');
  lines.push(`test.describe(${JSON.stringify(feature.name)}, () => {`);
  lines.push('  const session = feature(test);');
  lines.push(...body);
  lines.push('});');
  lines.push('');
  return {feature, outFile, code: lines.join('\n'), diagnostics};
}

type EmitStep = (step: StepModel, indent: string, state: ScenarioState) => string[];

function tagOptions(tags: string[]): string {
  const unique = [...new Set(tags)].filter((t) => t.startsWith('@'));
  return unique.length > 0 ? `, {tag: [${unique.map((t) => JSON.stringify(t)).join(', ')}]}` : '';
}

function emitScenario(scenario: ScenarioModel, feature: FeatureModel, uniqueTitle: (s: string) => string, emitStep: EmitStep): string[] {
  const out: string[] = [];
  const state: ScenarioState = {};
  out.push(`  test(${JSON.stringify(uniqueTitle(scenario.name))}${tagOptions([...feature.tags, ...scenario.tags])}, async ({browser}) => {`);
  out.push('    const page = await session.page(browser);');
  for (const step of [...feature.background, ...scenario.steps])
    out.push(...emitStep(step, '    ', state));
  out.push('  });');
  return out;
}

/** One test for the feature: the Background, then every scenario as a soft step on the same state. */
function emitJourney(feature: FeatureModel, uniqueTitle: (s: string) => string, emitStep: EmitStep, helpers: Set<string>): string[] {
  helpers.add('journey');
  const out: string[] = [];
  const state: ScenarioState = {};
  out.push(`  test(${JSON.stringify(feature.name)}${tagOptions([...feature.tags, ...feature.scenarios.flatMap((s) => s.tags)])}, async ({browser}) => {`);
  out.push('    const page = await session.page(browser);');
  out.push(`    const run = journey(test, ${feature.scenarios.length});`);
  for (const step of feature.background)
    out.push(...emitStep(step, '    ', state));
  for (const scenario of feature.scenarios) {
    out.push(`    await run.scenario(${JSON.stringify(uniqueTitle(scenario.name))}, async () => {`);
    for (const step of scenario.steps)
      out.push(...emitStep(step, '      ', state));
    out.push('    });');
  }
  out.push('    run.finish();');
  out.push('  });');
  return out;
}
