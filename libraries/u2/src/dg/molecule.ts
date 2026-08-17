/* Molecule support: bridged, never ported (OBJECTS.md). The sketcher is a plugin ecosystem
   behind the platform's `role: valueEditor` registration; u2 wraps its editor through the
   input bridge and its depiction through drawMolecule, and contains no chemistry itself. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Input} from '../core/input-base.js';
import {ObjectRenderer} from '../core/object-renderer.js';
import {fromDartInput} from './from-dart-input.js';
import {Editors} from './editors.js';

/** The platform's molecule editor (sketcher-backed, whatever backend the deployment runs) as a
 * u2 `Input<string>` over the molecule string — form validity, label alignment and the value
 * signal included. Requires the Chem package, like `ui.input.molecule` itself. */
export function moleculeInput(label: string, options?: {name?: string}): Input<string> {
  return fromDartInput<string>(ui.input.molecule(label), options?.name);
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
  create: (prop, options) => moleculeInput(options.label ?? prop.name, {name: options.name}),
});
