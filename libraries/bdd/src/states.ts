/* The states an element can be claimed to be in — one list for the `{state}` parameter type, the
   assertions and the VS Code settings `init` writes. */
export type State = 'visible' | 'hidden' | 'present' | 'absent' | 'enabled' | 'disabled' | 'checked' |
  'unchecked' | 'partially checked' | 'selected' | 'empty' | 'expanded' | 'collapsed' | 'focused' | 'invalid' | 'valid';

// "partially checked" before "checked": the alternation takes the first match
export const STATES: State[] = ['visible', 'hidden', 'present', 'absent', 'enabled', 'disabled', 'partially checked', 'checked',
  'unchecked', 'selected', 'empty', 'expanded', 'collapsed', 'focused', 'invalid', 'valid'];
