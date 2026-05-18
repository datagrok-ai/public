import * as DG from 'datagrok-api/dg';
import * as jsyaml from 'js-yaml';

import {_package} from '../package';

interface PromptSuggestion {
  label?: string;
  prompt: string;
}

type Scope = 'column' | 'viewer' | 'view' | 'global';

interface ScopeBlock {
  scope: Scope;
  key: string;
  match?: {semType?: string; type?: string};
  label: string;
  icon?: string;
  suggestions: PromptSuggestion[];
}

interface SuggestionsFile {
  powerSearch: Record<string, PromptSuggestion[]>;
  context: ScopeBlock[];
}

export interface ResolvedScope {
  key?: string;
  label: string;
  icon?: string;
  suggestions: PromptSuggestion[];
}

let _cache: SuggestionsFile | null = null;

async function loadSuggestions(): Promise<SuggestionsFile> {
  if (_cache)
    return _cache;
  const text = await _package.files.readAsText('suggestions.yaml');
  _cache = jsyaml.load(text) as SuggestionsFile;
  return _cache;
}

export async function resolveContextScopes(view: DG.ViewBase | null): Promise<ResolvedScope[]> {
  const tv = view instanceof DG.TableView ? view : null;
  const blocks = (await loadSuggestions()).context;
  const out: ResolvedScope[] = [];
  for (const b of blocks) {
    const hasCol = b.scope === 'column' && !!tv?.dataFrame?.columns.bySemType(b.match?.semType ?? '');
    const ok =
      b.scope === 'global' ||
      hasCol ||
      (b.scope === 'viewer' && !!tv && Array.from(tv.viewers).some((v) => v.type === b.match?.type)) ||
      (b.scope === 'view' && view?.type === b.match?.type);
    if (ok)
      out.push({key: b.key, label: b.label, icon: b.icon, suggestions: b.suggestions});
  }
  return out;
}

export async function loadPowerSearchScopes(): Promise<ResolvedScope[]> {
  const data = (await loadSuggestions()).powerSearch;
  return Object.entries(data).map(([label, suggestions]) => ({label, suggestions}));
}

export function showSuggestionsMenu(
  scopes: ResolvedScope[],
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
