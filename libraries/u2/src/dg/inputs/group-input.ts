/* Platform group and user picker over the plain TypeAhead: the same look-up the platform's Share
   dialog and `ui.input.userGroups` run (`groups/lookup` — groups, and users through their personal
   groups), a two-line row (name over the kind). The sibling of `userInput`. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {TypeAhead} from '../../components/inputs/typeahead.js';

export interface GroupInputOptions {
  placeholder?: string;
  /** Keeps only the groups it accepts. */
  accept?: (group: DG.Group) => boolean;
  /** Every pick, with the label it was shown under (a duplicate name carries the id's start),
   * after which the box is cleared for the next one — an "add" control. */
  onPick?: (group: DG.Group, label: string) => void;
}

export function groupInput(options: GroupInputOptions = {}): TypeAhead<DG.Group> {
  const accept = options.accept;
  const duplicates = new Set<string>();
  const label = (group: DG.Group): string =>
    duplicates.has(group.friendlyName) ? `${group.friendlyName} (${group.id.slice(0, 8)})` : group.friendlyName;
  const input = new TypeAhead<DG.Group>({
    source: async (query) => {
      const found = query === '' ? [] : await grok.dapi.groups.getGroupsLookup(query);
      const kept = accept === undefined ? found : found.filter(accept);
      duplicates.clear();
      const seen = new Set<string>();
      for (const g of kept) {
        if (seen.has(g.friendlyName))
          duplicates.add(g.friendlyName);
        seen.add(g.friendlyName);
      }
      return kept;
    },
    itemText: (group) => group.friendlyName,
    render: (group) => renderGroup(group, label(group)),
    placeholder: options.placeholder ?? 'Group or user…',
    minChars: 1,
  });
  const onPick = options.onPick;
  if (onPick !== undefined) {
    input.effect(() => {
      const group = input.selected.value;
      if (group === null)
        return;
      onPick(group, label(group));
      input.selected.value = null;
    });
  }
  return input;
}

function renderGroup(group: DG.Group, label: string): HTMLElement {
  const row = document.createElement('div');
  row.className = 'u2-typeahead-user';
  const text = document.createElement('div');
  text.className = 'u2-typeahead-user-text';
  const kind = group.personal ? 'user' : 'group';
  text.append(line('u2-typeahead-user-name', group.friendlyName),
    line('u2-typeahead-user-secondary', label === group.friendlyName ? kind : `${kind} · ${group.id.slice(0, 8)}`));
  row.append(text);
  return row;
}

function line(cls: string, text: string): HTMLElement {
  const el = document.createElement('div');
  el.className = cls;
  el.textContent = text;
  return el;
}
