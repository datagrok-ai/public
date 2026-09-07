/* The three vocabularies a feature file binds to — steps (verbs), elements and kinds (nouns),
   datasets (data) — plus the custom parameter types that carry them through cucumber
   expressions. Binding modules under `bindings/` call these at import time; the compiler and the
   runtime import the same modules, so both sides read one registry.
   Element names are global and unique: bindings/platform registers the Datagrok shell (toolbox,
   browse tab, grid …) and those names are reserved everywhere; an app registers its own names on a
   context — a named region whose vocabulary applies once a step `enters` it. */
import type {Page} from '@playwright/test';

export type StepKeyword = 'Given' | 'When' | 'Then' | 'Step';
export type StepFn = (page: Page, ...args: any[]) => Promise<void> | void;

export interface StepMeta {
  /** `ui` drives the real interface, `api` goes through `page.evaluate` — informational. */
  tier?: 'ui' | 'api';
  description?: string;
  /** The context (see `context()`) whose vocabulary applies to the steps after this one. */
  enters?: string;
  /** Back to the platform vocabulary after this step. */
  leaves?: boolean;
}

export interface StepDef {
  keyword: StepKeyword;
  expression: string;
  fn: StepFn;
  meta: StepMeta;
}

export const steps: StepDef[] = [];

function definer(keyword: StepKeyword) {
  return <F extends StepFn>(expression: string, fn: F, meta: StepMeta = {}): F => {
    steps.push({keyword, expression, fn, meta});
    return fn;
  };
}

/** The keyword is informational (Cucumber matches any keyword against any definition); spelled
 * this way so the Cucumber VS Code extension indexes the definitions. */
export const Given = definer('Given');
export const When = definer('When');
export const Then = definer('Then');
export const Step = definer('Step');

export interface ElementGestures {
  /** `mouse` clicks through the real pointer at the element's centre — canvases ignore DOM clicks. */
  click?: 'mouse' | 'dom';
  /** `keyboard` types key by key (Dart change listeners ignore `fill`). */
  type?: 'keyboard' | 'fill';
  /** What opens the element's popup: Dart comboboxes open on `mousedown`. */
  open?: 'mousedown' | 'click';
}

export interface ElementDef {
  /** Playwright selector, resolved inside the scope the phrase names (or the whole page). */
  selector: string;
  /** Narrows `selector` to elements whose text matches. */
  text?: string | RegExp;
  /** Scope used when the phrase names none: another registered element. */
  in?: string;
  aliases?: string[];
  gestures?: ElementGestures;
  /** Named sub-elements, addressable as `<part> of <element>`. */
  parts?: Record<string, string>;
  description?: string;
}

export interface ElementEntry extends ElementDef {
  name: string;
  /** The context the name belongs to — resolved inside its root first. Unset for global names. */
  context?: ContextEntry;
}

/** A named region with a vocabulary of its own. The context is a global element itself ("MSA
 * workbench should be visible" works anywhere); the names registered on it apply to the steps
 * after one declared with `enters`, and cannot repeat a global (platform) name. */
export interface ContextEntry extends ElementEntry {
  element(name: string, def: ElementDef): ElementEntry;
  alias(name: string, target: string): void;
  lookup(phrase: string): ElementEntry | undefined;
  elements(): ElementEntry[];
}

/** `exact-text` is Playwright's own text engine (the innermost element whose text is the qualifier),
 * independent of the kind's selector. */
export type MatchStrategy = 'name' | 'label' | 'text' | 'exact-text' | 'title' | 'aria' | 'placeholder' | 'dart';

export interface KindDef {
  /** Base locator for candidates of this kind. */
  selector: string;
  aliases?: string[];
  /** How the qualifier ("sequence column" in "sequence column input") narrows the candidates,
   * in priority order — the first strategy that finds something wins. */
  match: MatchStrategy[];
  /** Where the kind keeps its caption (the input label, the dialog title) for `label`/`title`. */
  labelSelector?: string;
  /** Dart `name=` conventions for the qualifier; `{q}` is the qualifier with spaces as dashes.
   * Compared case-insensitively. */
  dartNames?: string[];
  gestures?: ElementGestures;
  /** The editable control inside the element, for typing and value assertions. */
  editorSelector?: string;
  /** Named sub-elements of every element of the kind, addressable as `<part> of <element>`. */
  parts?: Record<string, string>;
  description?: string;
}

export interface KindEntry extends KindDef {
  name: string;
}

export interface DatasetDef {
  /** Platform location (`System:DemoFiles/demog.csv`). */
  path: string;
  aliases?: string[];
  description?: string;
}

export interface DatasetEntry extends DatasetDef {
  name: string;
}

export interface ParameterTypeDef {
  name: string;
  regexp: RegExp;
  transformer?: (text: string) => unknown;
  description?: string;
}

const elementByName = new Map<string, ElementEntry>();
const contextByName = new Map<string, ContextEntry>();
const kindByName = new Map<string, KindEntry>();
const datasetByName = new Map<string, DatasetEntry>();
export const parameterTypes: ParameterTypeDef[] = [];

export const ARTICLES = /^(the|a|an)\s+/;

/** Lowercase, single-spaced, article-free — the key both registration and lookup use. */
export function normalizePhrase(s: string): string {
  return s.trim().toLowerCase().replace(/\s+/g, ' ').replace(ARTICLES, '').replace(/[.,;:!]+$/, '');
}

function keysOf(name: string, aliases?: string[]): string[] {
  return [normalizePhrase(name), ...(aliases ?? []).map(normalizePhrase)];
}

function describeTaken(taken: ElementEntry): string {
  return taken.context ? `"${taken.name}" of ${taken.context.name}` : `element "${taken.name}" [${taken.selector}]`;
}

function register(entry: ElementEntry, keys: string[]): void {
  for (const key of keys) {
    const taken = elementByName.get(key);
    if (taken)
      throw new Error(`element "${key}" is already registered as ${describeTaken(taken)} — names are global; rename it, or register app names on a context`);
  }
  for (const key of keys)
    elementByName.set(key, entry);
}

export function element(name: string, def: ElementDef): ElementEntry {
  const entry: ElementEntry = {...def, name: normalizePhrase(name)};
  register(entry, keysOf(name, def.aliases));
  return entry;
}

/** Another name for a registered element; the target may be registered later. */
export function alias(name: string, target: string): void {
  const entry = elementByName.get(normalizePhrase(target));
  if (!entry)
    throw new Error(`alias "${name}": element "${target}" is not registered (register it first)`);
  register(entry, [normalizePhrase(name)]);
}

export function context(name: string, def: ElementDef): ContextEntry {
  const ctxName = normalizePhrase(name);
  const local = new Map<string, ElementEntry>();
  const registerLocal = (entry: ElementEntry, keys: string[]) => {
    for (const key of keys) {
      const taken = elementByName.get(key) ?? local.get(key);
      if (taken)
        throw new Error(`"${key}" is already registered as ${describeTaken(taken)} — platform names are reserved inside "${ctxName}" too; rename it`);
    }
    for (const key of keys)
      local.set(key, entry);
  };
  const ctx: ContextEntry = {
    ...def,
    name: ctxName,
    element(n, d) {
      const entry: ElementEntry = {...d, name: normalizePhrase(n), context: ctx};
      registerLocal(entry, keysOf(n, d.aliases));
      return entry;
    },
    alias(n, target) {
      const entry = local.get(normalizePhrase(target));
      if (!entry)
        throw new Error(`alias "${n}": "${target}" is not registered on ${ctxName} (register it first)`);
      registerLocal(entry, [normalizePhrase(n)]);
    },
    lookup: (phrase) => local.get(normalizePhrase(phrase)),
    elements: () => [...new Set(local.values())],
  };
  register(ctx, keysOf(name, def.aliases));
  for (const key of keysOf(name, def.aliases))
    contextByName.set(key, ctx);
  return ctx;
}

export function kind(name: string, def: KindDef): KindEntry {
  const entry: KindEntry = {...def, name: normalizePhrase(name)};
  for (const key of keysOf(name, def.aliases))
    kindByName.set(key, entry);
  return entry;
}

export function dataset(name: string, def: DatasetDef): DatasetEntry {
  const entry: DatasetEntry = {...def, name: normalizePhrase(name)};
  for (const key of keysOf(name, def.aliases))
    datasetByName.set(key, entry);
  return entry;
}

export function defineParameterType(def: ParameterTypeDef): ParameterTypeDef {
  parameterTypes.push(def);
  return def;
}

export function lookupElement(phrase: string): ElementEntry | undefined {
  return elementByName.get(normalizePhrase(phrase));
}

export function lookupContext(phrase: string): ContextEntry | undefined {
  return contextByName.get(normalizePhrase(phrase));
}

export function lookupKind(phrase: string): KindEntry | undefined {
  return kindByName.get(normalizePhrase(phrase));
}

export function lookupDataset(name: string): DatasetEntry | undefined {
  return datasetByName.get(normalizePhrase(name));
}

/** Every registered element name and alias, longest first. */
export function elementNames(): string[] {
  return [...elementByName.keys()].sort((a, b) => b.length - a.length || a.localeCompare(b));
}

/** Every registered kind name and alias, longest first — suffix matching needs that order. */
export function kindNames(): string[] {
  return [...kindByName.keys()].sort((a, b) => b.length - a.length || a.localeCompare(b));
}

export function datasetNames(): string[] {
  return [...datasetByName.keys()].sort();
}

/** Unique registered elements/contexts/kinds/datasets (aliases collapsed). */
export function allElements(): ElementEntry[] {
  return [...new Set(elementByName.values())];
}

export function allContexts(): ContextEntry[] {
  return [...new Set(contextByName.values())];
}

export function allKinds(): KindEntry[] {
  return [...new Set(kindByName.values())];
}

export function allDatasets(): DatasetEntry[] {
  return [...new Set(datasetByName.values())];
}

export function resetRegistry(): void {
  steps.length = 0;
  parameterTypes.length = 0;
  elementByName.clear();
  contextByName.clear();
  kindByName.clear();
  datasetByName.clear();
}
