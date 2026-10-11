/** Short notices shown wherever the names of the original programs appear:
 * this package reimplements HMMER 3.4 (translated to Rust, compiled to
 * WebAssembly) and ANARCI (translated to TypeScript); neither original
 * program runs here, and the package is not affiliated with their authors. */
import * as ui from 'datagrok-api/ui';

export const HMMER_NOTE = 'HMMER 3.4 reimplemented in Rust and compiled to WebAssembly (Rusty-HMMER), ' +
  'not the original HMMER; results are identical.';
export const ANARCI_NOTE = 'ANARCI reimplemented in TypeScript on Rusty-HMMER (HMMER 3.4 in Rust and ' +
  'WebAssembly), not the original programs; results are identical.';
export const CLAN_NOTE = 'Clan competition is this package\'s implementation of Pfam\'s rule, checked against ' +
  'Pfam\'s pfam_scan.pl.';

/** A small muted note for dialogs and panels. */
export function noteDiv(...notes: string[]): HTMLElement {
  const div = ui.divText(notes.join(' '), {style: {fontSize: '11px', lineHeight: '1.35', color: 'var(--grey-5)',
    maxWidth: '440px', marginTop: '10px'}});
  div.classList.add('hmm-reimplementation-note');
  return div;
}
