/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';

import $ from 'cash-dom';

import {extractAtomDataV3000} from './mol-transformations';

/** Draw molecule on the canvas and append it to the specified div, with the
 * option of zoom-in */
export async function drawMolecule(
  canvasWidth: number, canvasHeight: number,
  moleculeImgDiv: HTMLDivElement,
  molfile: string
): Promise<void> {
  // clear the div's content if any
  moleculeImgDiv.innerHTML = '';

  const canvas = ui.canvas(canvasWidth * window.devicePixelRatio, canvasHeight * window.devicePixelRatio);

  async function drawMolfileOnCanvas(canvas: HTMLCanvasElement): Promise<void> {
    await grok.functions.call('Chem:canvasMol', {
      x: 0, y: 0, w: canvas.width, h: canvas.height, canvas: canvas,
      molString: molfile, scaffoldMolString: '',
      options: {normalizeDepiction: false, straightenDepiction: false}
    });
  }

  async function drawZoomedInMolecule(): Promise<void> {
    try {
      const addDiv = ui.div([], {style: {overflowX: 'scroll'}});

      // addDiv size required, but now available before dialog show()
      const atomCoordinates = extractAtomDataV3000(molfile);
      // const cw: number = $(window).width() * 0.80; // addDiv.clientWidth
      const clientHeight: number = $(window).height() * 0.70; // addDiv.clientHeight
      const molWidth: number = Math.max(...atomCoordinates.x) - Math.min(...atomCoordinates.x);
      const molHeight: number = Math.max(...atomCoordinates.y) - Math.min(...atomCoordinates.y);

      // const wR: number = cw / molWidth;
      const hR: number = clientHeight / molHeight;
      const r: number = hR; // Math.max(wR, hR);
      const dialogCanvasWidth = r * molWidth;
      const dialogCanvasHeight = r * molHeight;

      const dialogCanvas = ui.canvas(
        dialogCanvasWidth * window.devicePixelRatio, dialogCanvasHeight * window.devicePixelRatio
      );
      dialogCanvas.style.width = `${dialogCanvasWidth}px`;
      dialogCanvas.style.height = `${dialogCanvasHeight}px`;
      await drawMolfileOnCanvas(dialogCanvas);

      addDiv.appendChild(dialogCanvas);
      ui.dialog('Molecule')
        .add(addDiv)
        .showModal(true);
    } catch (err) {
      const errStr = errorToConsole(err);
      console.error(errStr);
    }
  };

  // Draw zoomed-out molecule
  canvas.style.width = `${canvasWidth}px`;
  canvas.style.height = `${canvasHeight}px`;
  canvas.style.borderStyle = 'solid';
  canvas.style.borderColor = 'blue';
  drawMolfileOnCanvas(canvas);

  // Dialog with zoomed-in molecule
  // canvas.addEventListener('click', drawZoomedInMolecule);
  $(canvas).on('click', drawZoomedInMolecule);
  $(canvas).on('mouseover', () => $(canvas).css('cursor', 'grab')); // for some reason 'zoom-in' value wouldn't work
  $(canvas).on('mouseout', () => $(canvas).css('cursor', 'default'));

  moleculeImgDiv.append(canvas);
}
