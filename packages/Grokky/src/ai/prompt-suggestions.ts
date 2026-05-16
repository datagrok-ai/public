import * as DG from 'datagrok-api/dg';
import * as jsyaml from 'js-yaml';

import {_package} from '../package';

interface PromptSuggestion {
  label?: string;
  prompt: string;
}

interface SuggestionsFile {
  powerSearch: Record<string, PromptSuggestion[]>;
  context: Record<string, PromptSuggestion[]>;
}

let _cache: SuggestionsFile | null = null;

async function loadSuggestions(): Promise<SuggestionsFile> {
  if (_cache)
    return _cache;
  const text = await _package.files.readAsText('suggestions.yaml');
  _cache = jsyaml.load(text) as SuggestionsFile;
  return _cache;
}

export async function showSuggestionsMenu(
  section: keyof SuggestionsFile,
  onSelect: (prompt: string) => void,
) {
  const menu = DG.Menu.popup();
  for (const [category, items] of Object.entries((await loadSuggestions())[section])) {
    const group = menu.group(category);
    for (const s of items)
      group.item(s.label ?? s.prompt, () => onSelect(s.prompt));
  }
  menu.show();
}
