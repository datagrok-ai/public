/* Molecule support: bridged, never ported (OBJECTS.md). The sketcher is a plugin ecosystem
   behind the platform's `role: valueEditor` registration; u2 wraps its editor through the
   input bridge and its depiction through drawMolecule, and contains no chemistry itself. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Input, labelText} from '../../core/input-base.js';
import type {InputOptions} from '../../core/input-base.js';
import {ObjectRenderer} from '../../core/object-renderer.js';
import {fromDartInput} from './from-dart-input.js';
import {Editors} from '../forms/editors.js';

/** The platform's molecule editor (sketcher-backed, whatever backend the deployment runs) as a
 * u2 `Input<string>` over the molecule string — form validity, label alignment and the value
 * signal included. Requires the Chem package, like `ui.input.molecule` itself. */
export function moleculeInput(label: string,
  options: Pick<InputOptions<string>, 'name' | 'value' | 'onChanged' | 'inline' | 'nullable'> = {}): Input<string> {
  const {name, value, ...rest} = options;
  const dgInput = ui.input.molecule(label, {value: value ?? ''});
  const input = fromDartInput<string>(dgInput, name, rest);
  sketchInDialog(input, dgInput.root);
  return input;
}

/** What a structure the sketcher drew nothing for is told. */
const UNREADABLE_STRUCTURE = 'Cannot read the structure';

/** The platform input's dialog keeps what its sketcher drew and nothing else: a SMILES it could not
 * read vanished without a word. The same sketcher opens here instead, and text it drew nothing for
 * stays the value — invalid unless the platform reads it as SMILES — for the user to see and fix. */
function sketchInDialog(input: Input<string>, root: HTMLElement): void {
  const host = root.querySelector<HTMLElement>('.d4-input-molecule-canvas-host');
  if (host === null)
    return;
  let unreadable: string | null = null;
  input.addValidator((value) => value !== '' && value === unreadable ? UNREADABLE_STRUCTURE : null);
  // ahead of the platform's own click handler, which never runs
  const open = (e: Event) => {
    if (input.readOnly)
      return;
    e.stopImmediatePropagation();
    const current = input.value.peek() ?? '';
    const sketcher = new grok.chem.Sketcher();
    if (current === unreadable)
      sketcher.molInput.value = current;
    else
      sketcher.setMolecule(current);
    ui.dialog().add(sketcher.root).onOK(() => {
      const typed = sketcher.molInput.value.trim();
      const drawn = sketcher.isEmpty() ? '' : sketcher.getSmiles();
      unreadable = drawn === '' && typed !== '' && !grok.chem.isMolBlock(typed) && !grok.chem.checkSmiles(typed) ?
        typed : null;
      input.value.value = drawn !== '' ? drawn : typed;
      input.touch();
      input.revalidate();
    }).show();
  };
  host.addEventListener('click', open, {capture: true});
  input.own(() => host.removeEventListener('click', open, {capture: true}));
}

export interface MoleculeRendererOptions {
  width?: number;
  height?: number;
}

/** Renders molecule strings (SMILES / molblock) on collection surfaces — typeahead rows, chips,
 * cards. `drawMolecule` returns its host synchronously and fills it when Chem delivers, which
 * is exactly the placeholder-swap contract collection controls expect. */
export function moleculeRenderer(options: MoleculeRendererOptions = {}): ObjectRenderer<string> {
  const w = options.width ?? 150;
  const h = options.height ?? 100;
  const draw = (s: string, dw: number, dh: number) => {
    const host = grok.chem.drawMolecule(s, dw, dh);
    host.style.width = `${dw}px`;
    host.style.height = `${dh}px`;
    return host;
  };
  return {
    caption: (s) => s,
    listItem: (s) => draw(s, w, h),
    markup: (s) => draw(s, w, h),
    tooltip: (s) => draw(s, Math.max(w, 250), Math.max(h, 160)),
  };
}

Editors.register({
  match: (prop) => prop.semType === 'Molecule',
  create: (prop, options) => moleculeInput(options.inline ? '' : labelText(options.label) ?? prop.name!,
    {name: options.name, value: options.value, onChanged: options.onChanged, inline: options.inline,
      nullable: options.nullable}),
});
