/* The `molecules` tier: Crux Sketch, Chem's own molecule sketcher (Chem's src/crux/crux-sketcher.ts), as features of
   any package drive it wherever a host opens it (a cell editor, a molecule input's dialog, an inline host in an app).
   Its root is named `cruxSketch` (data-u2-name). The sketcher reports its status to the platform (getWidgetStatus): a
   hit area for every atom and bond it draws, "atom 0", "bond 0"; its parts, "canvas", "actions" (the top toolbar),
   "tools" (the left one), "elements" (the element palette), "templates" (the ring bar) and "label editor" while it is
   open; "tool <name>" for each control a toolbar shows, by Crux's name for it ("tool ring.benzene", "tool undo"); and
   the readings "ready", "pending", "smiles" (Crux's own), "atoms", "bonds", "mode", "selected atoms", "selected bonds",
   "tool", "query", "empty" and "changes". The viewers tier's area and reading steps take it as "crux sketcher widget",
   and its gestures settle on its isRenderPending and onRendered. The same parts are elements here (`<part> of crux
   sketcher widget`), and Crux's own controls are found by their test ids (crux-sketch's docs/conventions/test-ids.md),
   in its open shadow root. A feature that draws on them pins Crux (`the molecule sketcher is "Crux"`) and is tagged
   @sketcher-controls. */
import {element} from '../../../src/registry.js';

const ROOT = '[data-u2-name="cruxSketch"]';

element('crux sketcher widget', {selector: ROOT,
  parts: {
    'canvas': 'crux-sketch [data-testid="canvas"]',
    'actions': 'crux-sketch [data-testid="toolbar.top"]',
    'tools': 'crux-sketch [data-testid="toolbar.left"]',
    'elements': 'crux-sketch [data-testid="toolbar.right"]',
    'templates': 'crux-sketch [data-testid="toolbar.bottom"]',
    'label editor': 'crux-sketch [data-testid="canvas.label-editor"]',
  },
  description: 'Crux wherever a host shows it (data-u2-name cruxSketch), its status the platform\'s: areas "atom N", ' +
    '"bond N", canvas, actions, tools, elements, templates, label editor and "tool <name>"; readings ready, pending, ' +
    'smiles, atoms, bonds, mode, selected atoms, selected bonds, tool, query, empty, changes'});

/** Crux's own controls by their test ids. */
const PARTS: [string, string, string][] = [
  ['crux canvas', 'canvas', 'the drawing area of Crux'],
  ['crux preview', 'overlay.preview', 'what Crux shows under the pointer before a click places it (a ring\'s preview)'],
  ['crux single bond tool', 'toolbar.bond.single', 'the single bond tool of Crux\'s left toolbar'],
  ['crux benzene tool', 'toolbar.ring.benzene', 'the benzene ring of Crux\'s ring bar'],
  ['crux select tool', 'toolbar.select', 'the selection palette\'s button (the rectangle at first)'],
  ['crux clear button', 'toolbar.clear', 'Clear canvas, on Crux\'s top toolbar'],
  ['crux undo button', 'toolbar.undo', 'Undo, on Crux\'s top toolbar'],
  ['crux flip horizontal button', 'canvas.selection-actions.flip-horizontal', 'the flip button beside a selection'],
  ['crux nitrogen tool', 'toolbar.element.n', 'the nitrogen of Crux\'s element palette: a click on an atom makes it N'],
  ['crux R-group tool', 'toolbar.rgroup', 'the R-group palette\'s button, the R-group label tool at first: a click on an atom opens the R-Group dialog'],
  ['crux R1 button', 'dialog.rgroup.r.1', 'R1 in Crux\'s R-Group dialog'],
  ['crux R2 button', 'dialog.rgroup.r.2', 'R2 in Crux\'s R-Group dialog'],
  ['crux R-Group OK button', 'dialog.rgroup.ok', 'OK of Crux\'s R-Group dialog'],
];
for (const [name, id, description] of PARTS)
  element(name, {selector: `${ROOT} crux-sketch [data-testid="${id}"]`, description});
