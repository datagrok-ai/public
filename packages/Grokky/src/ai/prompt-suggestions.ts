import * as DG from 'datagrok-api/dg';
import * as jsyaml from 'js-yaml';

import {_package} from '../package';

type Scope = 'column' | 'view' | 'global';
type Context = 'powerSearch' | 'panel';

export interface Suggestion {
  label?: string;
  prompt: string;
  prefer?: Record<string, string>;
}

export interface Block {
  context: Context;
  scope: Scope;
  key: string;
  match?: {semType?: string; type?: string};
  label: string;
  icon?: string;
  suggestions: Suggestion[];
}

// TODO: LLM-curated suggestions vs current human + view-agnostic heuristics?
const SLOT = /\{(\w+):([^}]+)}/g;
let _cache: Block[] | null = null;

async function load(): Promise<Block[]> {
  if (!_cache)
    _cache = jsyaml.load(await _package.files.readAsText('suggestions.yaml')) as Block[];
  return _cache;
}

function colMatches(c: DG.Column, filter: string): boolean {
  if (filter === 'numerical')
    return c.isNumerical;
  if (filter === 'categorical')
    return c.isCategorical;
  if (filter.startsWith('semType='))
    return c.semType === filter.slice(8);
  return false;
}

function resolveSuggestion(s: Suggestion, df: DG.DataFrame | null): Suggestion | null {
  const picked: Record<string, string> = {};
  const taken = new Set<string>();
  let cols: DG.Column[] | null = null;

  for (const [, name, filter] of `${s.label ?? ''}\n${s.prompt}`.matchAll(SLOT)) {
    if (name in picked) continue;
    if (!df) return null;
    if (!cols) cols = df.columns.toList();
    const candidates = cols.filter((c) => !taken.has(c.name) && colMatches(c, filter));
    if (!candidates.length) return null;

    const prefer = s.prefer?.[name];
    const re = prefer ? new RegExp(prefer.match(/name~="([^"]+)"/)?.[1] ?? prefer, 'i') : null;
    const pick = (re && candidates.find((c) => re.test(c.name))) || candidates[0];
    picked[name] = pick.name;
    taken.add(pick.name);
  }

  const fill = (t: string): string => t.replace(SLOT, (_, n) => picked[n]);
  return {label: s.label ? fill(s.label) : undefined, prompt: fill(s.prompt)};
}

export async function resolveScopes(context: Context, view: DG.ViewBase | null = null): Promise<Block[]> {
  const df = view instanceof DG.TableView ? view.dataFrame ?? null : null;
  const out: Block[] = [];
  for (const b of await load()) {
    if (b.context !== context) continue;
    const applies = b.scope === 'global' ||
      (b.scope === 'view' && view?.type === b.match?.type) ||
      (b.scope === 'column' && !!df?.columns.bySemType(b.match?.semType ?? ''));
    if (!applies) continue;

    const suggestions = b.suggestions
      .map((s) => resolveSuggestion(s, df))
      .filter((s): s is Suggestion => !!s);
    if (suggestions.length)
      out.push({...b, suggestions});
  }
  return out;
}

export function showSuggestionsMenu(
  scopes: Block[],
  onSelect: (prompt: string) => void,
  causedBy?: MouseEvent,
): void {
  const menu = DG.Menu.popup();
  for (const s of scopes) {
    const group = menu.group(s.label);
    for (const sg of s.suggestions)
      group.item(sg.label ?? sg.prompt, () => onSelect(sg.prompt));
  }
  menu.show(causedBy ? {causedBy} : undefined);
}
