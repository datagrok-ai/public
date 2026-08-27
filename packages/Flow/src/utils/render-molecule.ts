/** Molecule → widget for the Flow catalog. The molecule is drawn ONCE on a large
 *  square backing canvas and CSS scales it down, so the in-node preview (whose
 *  container the user resizes freely) never shows an upscaled, blurry bitmap. */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {setTid} from './test-ids';

/** Backing-store resolution (square, aspect ratio 1). Well above the 900px
 *  in-node preview cap, so the drawing outresolves any display size. */
export const MOLECULE_CANVAS_SIZE = 1000;

export async function renderMoleculeWidget(molecule: string): Promise<DG.Widget> {
  // Not ui.canvas — it pins the CSS size to the pixel size, and here they must differ.
  const canvas = document.createElement('canvas');
  canvas.width = MOLECULE_CANVAS_SIZE;
  canvas.height = MOLECULE_CANVAS_SIZE;
  setTid(canvas, 'render-molecule-canvas');
  if (molecule?.trim())
    await grok.chem.canvasMol(0, 0, MOLECULE_CANVAS_SIZE, MOLECULE_CANVAS_SIZE, canvas, molecule);
  // AFTER the draw: Chem pins style.width/height to the logical size (and scales the
  // backing store by devicePixelRatio). We scale by CSS instead. Inline, not a
  // funcflow.css class — the widget renders wherever the function is called
  // (dashboards, other packages), where Flow's stylesheet may not be loaded.
  canvas.style.display = 'block';
  canvas.style.width = '100%';
  canvas.style.height = 'auto';
  canvas.style.aspectRatio = '1';
  const root = ui.div([canvas], 'ff-render-molecule');
  setTid(root, 'render-molecule');
  return DG.Widget.fromRoot(root);
}
