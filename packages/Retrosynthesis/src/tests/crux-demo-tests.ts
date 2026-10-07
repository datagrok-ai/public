import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';
import {DEMO_MOLECULE} from '../const';

/* The demo's pane recomputes on every change of its sketcher (HOST-071): with Crux, Chem's own sketcher, as the
   session's, a completed edit asks the pane for one recomputation, and opening it asks once for the demo's molecule.
   The pane's function is stubbed: the recomputation itself (the AiZynthFinder container behind it) is not what is
   claimed, and the test needs no container. */

interface CruxElement extends HTMLElement {
  readonly smiles: string;
  readonly isPending: boolean;
  /** Where Crux draws each atom and bond (crux-sketch API-056), in px from the element's corner. */
  readonly positions: {atoms: ({x: number, y: number} | null)[], bonds: ({x: number, y: number} | null)[]};
}

/** Long enough for the demo's handler and for any further event a change could still cause. */
const SETTLE_MS = 300;

/** A press and release on Crux's canvas at a point of the element, as the pointer's events reach the canvas. */
function pressAt(crux: CruxElement, p: {x: number, y: number}): void {
  const canvas = crux.shadowRoot!.querySelector('[data-testid="canvas"]')!;
  const r = crux.getBoundingClientRect();
  const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true,
    clientX: r.left + p.x, clientY: r.top + p.y, pointerId: 1, pointerType: 'mouse', isPrimary: true, button: 0,
    buttons});
  canvas.dispatchEvent(new PointerEvent('pointerdown', init(1)));
  canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
}

/** Presses a button of Crux's toolbar, as the user does. */
function press(crux: CruxElement, testId: string): void {
  const button = crux.shadowRoot!.querySelector<HTMLElement>(`[data-testid="${testId}"]`);
  if (!button)
    throw new Error(`Crux has no ${testId}`);
  button.click();
}

category('Crux', () => {
  test('the demo asks its pane for one recomputation per completed edit in Crux', async () => {
    const sketcherWas = grok.chem.currentSketcherType;
    const paneWas = PackageFunctions.retroSynthesisPath;
    const asked: string[] = [];
    PackageFunctions.retroSynthesisPath = (molecule: string): DG.Widget => {
      asked.push(molecule);
      return new DG.Widget(ui.div());
    };
    grok.chem.currentSketcherType = 'Crux';
    try {
      await PackageFunctions.retrosynthesisDemo();
      const view = grok.shell.v;
      expect(view.name, 'Retrosynthesis Demo');
      let crux: CruxElement | null = null;
      await awaitCheck(() => {
        crux = view.root.querySelector<CruxElement>('.crux-sketcher crux-sketch');
        return crux !== null && crux.positions.atoms.length > 0 && !crux.isPending;
      }, 'Crux did not show the demo\'s molecule', 30000);
      // opening it: the molecule the demo sets, one change, the pane asked once (for the demo's own molecule)
      await awaitCheck(() => asked.length === 1, `the pane was asked ${asked.length} times on opening`, 5000);
      expect(asked[0], DEMO_MOLECULE);
      await delay(SETTLE_MS);
      expect(asked.length, 1);

      // a completed edit: a single bond drawn from atom 0
      const atoms = crux!.positions.atoms.length;
      press(crux!, 'toolbar.bond.single');
      pressAt(crux!, crux!.positions.atoms[0]!);
      await awaitCheck(() => crux!.positions.atoms.length === atoms + 1, 'the bond was not drawn', 5000);
      await awaitCheck(() => asked.length === 2, `after one edit the pane was asked ${asked.length} times`, 5000);
      await delay(SETTLE_MS);
      expect(asked.length, 2);
      expect(asked[1], crux!.smiles);

      // and a second: one more
      pressAt(crux!, crux!.positions.atoms[1]!);
      await awaitCheck(() => crux!.positions.atoms.length === atoms + 2, 'the second bond was not drawn', 5000);
      await awaitCheck(() => asked.length === 3, `after two edits the pane was asked ${asked.length} times`, 5000);
      await delay(SETTLE_MS);
      expect(asked.length, 3);
      expect(asked[2], crux!.smiles);
      view.close();
    } finally {
      PackageFunctions.retroSynthesisPath = paneWas;
      grok.chem.currentSketcherType = sketcherWas;
    }
  }, {timeout: 60000});
});
