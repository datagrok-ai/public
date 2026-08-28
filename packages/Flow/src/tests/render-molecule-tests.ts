/** Flow:renderMolecule — a widget with the molecule drawn once on a large square
 *  backing canvas (CSS scales it down), so the in-node preview stays crisp. */
import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {registerBuiltinNodes, registerAllFunctions, ensureFuncNodeType, createNode} from '../rete/node-factory';
import {FlowNode, supportsInlinePreview} from '../rete/scheme';
import {renderMoleculeWidget, MOLECULE_CANVAS_SIZE} from '../utils/render-molecule';

const BENZENE = 'c1ccccc1';

function canvasOf(w: DG.Widget): HTMLCanvasElement {
  const canvas = w.root.querySelector<HTMLCanvasElement>('[data-testid="ff-render-molecule-canvas"]');
  if (!canvas) throw new Error('no canvas in the widget');
  return canvas;
}

/** Distinct RGBA values over a coarse sample grid — a drawn molecule has several. */
function distinctColors(canvas: HTMLCanvasElement): number {
  const ctx = canvas.getContext('2d')!;
  const img = ctx.getImageData(0, 0, canvas.width, canvas.height).data;
  const seen = new Set<number>();
  const step = 16 * 4;
  for (let i = 0; i < img.length; i += step)
    seen.add((img[i] << 24) | (img[i + 1] << 16) | (img[i + 2] << 8) | img[i + 3]);
  return seen.size;
}

category('Flow: render molecule', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
  });

  test('registered, catalog-included, and running through the platform', async () => {
    const func = DG.Func.find({package: 'Flow', name: 'renderMolecule'})[0];
    expect(func != null, true, 'Flow:renderMolecule is registered');
    expect(String((func.options as {[k: string]: string})['includeInFlow']), 'true', 'opted into the catalog');
    const widget = await grok.functions.call('Flow:renderMolecule', {molecule: BENZENE}) as DG.Widget;
    expect(widget?.root != null, true, 'returns a widget with a root');
    const canvas = canvasOf(widget);
    // Chem scales the backing store by devicePixelRatio — at least the requested size.
    expect(canvas.width >= MOLECULE_CANVAS_SIZE, true, `backing width ${canvas.width}`);
    expect(canvas.height, canvas.width, 'square backing store (aspect ratio 1)');
    expect(distinctColors(canvas) > 1, true, 'the molecule was actually drawn');
  });

  test('CSS scales the big canvas down — crisp, square at any display size', async () => {
    const widget = await renderMoleculeWidget(BENZENE);
    const host = ui.div([widget.root], {style: {width: '300px', position: 'absolute', left: '-10000px'}});
    document.body.appendChild(host);
    try {
      const canvas = canvasOf(widget);
      expect(canvas.clientWidth, 300, 'display width follows the container');
      expect(canvas.clientHeight, 300, 'display height keeps aspect ratio 1');
      expect(canvas.width >= MOLECULE_CANVAS_SIZE, true, 'backing store stays full-resolution');
      expect(canvas.width > canvas.clientWidth, true, 'more pixels than displayed — never upscaled/blurry');
    } finally {
      host.remove();
    }
  });

  test('an empty molecule yields a blank widget, no throw', async () => {
    const widget = await renderMoleculeWidget('');
    const canvas = canvasOf(widget);
    expect(canvas.width, MOLECULE_CANVAS_SIZE);
    expect(distinctColors(canvas), 1, 'nothing drawn');
  });

  test('the node supports the in-node preview and takes a Molecule string', async () => {
    const func = DG.Func.find({package: 'Flow', name: 'renderMolecule'})[0];
    const node = createNode(ensureFuncNodeType(func)) as FlowNode;
    expect(node != null, true, 'node created');
    expect(supportsInlinePreview(node), true, 'widget output → in-node preview available');
    const molInput = node.dgFunc?.inputs.find((p) => p.name === 'molecule');
    expect(molInput?.semType ?? '', 'Molecule', 'input is semType Molecule (sketcher-compatible)');
    expect(molInput?.propertyType, 'string', 'as a string param');
  });
});
