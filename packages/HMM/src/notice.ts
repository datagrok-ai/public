/** Short notices, in tooltips wherever the names of the original programs
 * appear: this package reimplements HMMER 3.4 (translated to Rust, compiled
 * to WebAssembly) and ANARCI (translated to TypeScript); neither original
 * program runs here, and the package is not affiliated with their authors. */
import * as ui from 'datagrok-api/ui';

export const HMMER_NOTE = 'HMMER 3.4 reimplemented in Rust and compiled to WebAssembly (Rusty-HMMER), ' +
  'not the original HMMER; results are identical.';
export const ANARCI_NOTE = 'ANARCI reimplemented in TypeScript on Rusty-HMMER (HMMER 3.4 in Rust and ' +
  'WebAssembly), not the original programs; results are identical.';
export const CLAN_NOTE = 'Clan competition is this package\'s implementation of Pfam\'s rule, checked against ' +
  'Pfam\'s pfam_scan.pl.';

/** A tooltip text followed by notices (dialog inputs). */
export function withNotes(text: string, ...notes: string[]): string {
  return [text, ...notes].filter((s) => s).join(' ');
}

/** Show notices in the tooltip of `element` (panel content); returns it. */
export function noteTooltip<T extends HTMLElement>(element: T, ...notes: string[]): T {
  ui.tooltip.bind(element, notes.join(' '));
  element.classList.add('hmm-reimplementation-note');
  return element;
}
