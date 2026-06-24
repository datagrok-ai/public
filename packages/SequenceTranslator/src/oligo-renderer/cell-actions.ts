/**
 * Cell-level user actions for OligoNucleotide cells.
 *
 * The corresponding decorated entry points live in `package.ts` so the
 * platform can discover them — those methods are kept thin and forward to
 * the implementations below. Keeping the logic here means `package.ts` stays
 * structural and the oligo subsystem stays self-contained.
 *
 *   - `openOligoCanvasDialog`       — double-click handler: full-screen modal
 *     with a hi-res canvas rendering, hover interactions reused from the
 *     grid cell renderer.
 *   - `openOligoHelmEditorDialog`   — "Edit HELM" action: HELM Web Editor;
 *     OK writes the edited HELM back to the cell.
 *   - `copyHelmToClipboard`         — "Copy as HELM" action.
 *   - `copyDuplexImageToClipboard`  — "Copy as Image" action: hi-res PNG with
 *     a transparent background, trimmed to the duplex's natural extent.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {parseHelmDuplex} from './helm-parser';
import {drawDuplex, hitTest, DuplexLayout} from './canvas-renderer';
import {showMonomerTooltip} from './tooltip';

/** Reference layout dimensions used as the seed for sizing both the dialog
 * canvas and the clipboard image. Wide enough that any real-world duplex
 * fits without wrapping; the right edge is trimmed back to the real extent. */
const REF_W = 1500;
const REF_H = 120;
const RIGHT_PAD = 20;

/* -------------------------------------------------------------------------- *
 * Canvas viewer dialog (double-click on a cell)
 * -------------------------------------------------------------------------- */

export function openOligoCanvasDialog(cell: DG.GridCell): void {
  const helm = String(cell.cell?.value ?? '');
  if (!helm) return;
  const model = parseHelmDuplex(helm);

  // Probe layout at the reference dimensions so we discover the natural
  // horizontal extent of this particular duplex (varies with conjugates,
  // multi-char bases, strand length).
  const refLayout = drawDuplex(
    null as unknown as CanvasRenderingContext2D, 0, 0, REF_W, REF_H, model, undefined, true);
  const trimmedW = Math.max(REF_H, computeMaxChipEndX(refLayout) + RIGHT_PAD);

  // Scale the trimmed layout up to fill the modal — screen width minus 200px.
  const targetW = Math.max(400, window.innerWidth - 200);
  const scale = targetW / trimmedW;
  const finalW = trimmedW * scale; // === targetW by construction
  const finalH = REF_H * scale;

  // Backing canvas: dimensions in CSS px, pixel-density bumped by devicePixelRatio.
  const dpr = window.devicePixelRatio || 1;
  const canvas = document.createElement('canvas');
  canvas.width = Math.max(1, Math.round(finalW * dpr));
  canvas.height = Math.max(1, Math.round(finalH * dpr));
  canvas.style.width = `${finalW}px`;
  canvas.style.height = `${finalH}px`;
  canvas.style.display = 'block';
  const g = canvas.getContext('2d');
  if (!g) return;
  g.scale(dpr, dpr);
  const finalLayout = drawDuplex(g, 0, 0, finalW, finalH, model);

  // Hover support — same hit-test → tooltip pipeline the grid renderer uses,
  // just with canvas-local coords instead of GridCell.bounds-relative ones.
  canvas.addEventListener('mousemove', (e) => {
    const rect = canvas.getBoundingClientRect();
    const localX = e.clientX - rect.left;
    const localY = e.clientY - rect.top;
    const hit = hitTest(localX, localY, model, finalLayout);
    if (!hit) {
      ui.tooltip.hide();
      return;
    }
    showMonomerTooltip(hit, e.clientX, e.clientY);
  });
  canvas.addEventListener('mouseleave', () => ui.tooltip.hide());

  const container = ui.div([canvas], {style: {
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    width: '100%',
    height: '100%',
  }});

  ui.dialog({title: 'Oligonucleotide', showHeader: true, showFooter: false})
    .add(container)
    .show({modal: true, fullScreen: true});
}

/* -------------------------------------------------------------------------- *
 * HELM Web Editor action ("Edit HELM")
 * -------------------------------------------------------------------------- */

export async function openOligoHelmEditorDialog(value: DG.SemanticValue): Promise<void> {
  if (!value?.cell) return;
  const cell = value.cell;
  // Helm:editMoleculeCell can't be reused: it calls seqHelper.getSeqHandler(col),
  // which throws unless semType === Macromolecule. OligoNucleotide columns have
  // semType=OligoNucleotide, so we drive the HELM Web Editor directly.
  const helmHelper = await getHelmHelper();
  const view = ui.div();
  view.style.height = '100%';
  const app = helmHelper.createWebEditorApp(view, String(cell.value ?? ''));
  ui.dialog({showHeader: false, showFooter: true})
    .add(view)
    .onOK(() => {
      const helmValue = app.canvas!.getHelm(true)
        .replace(/<\/span>/g, '')
        .replace(/<span style='background:#bbf;'>/g, '');
      cell.value = helmValue;
    })
    .show({modal: true, fullScreen: true});
}

/* -------------------------------------------------------------------------- *
 * Copy as HELM
 * -------------------------------------------------------------------------- */

export function copyHelmToClipboard(value: DG.SemanticValue): void {
  const helm = String(value?.value ?? '');
  if (!helm) return;
  navigator.clipboard.writeText(helm).then(() => grok.shell.info('HELM copied to clipboard'));
}

/* -------------------------------------------------------------------------- *
 * Copy as Image
 * -------------------------------------------------------------------------- */

const IMAGE_SCALE = 8;

export function copyDuplexImageToClipboard(value: DG.SemanticValue): void {
  if (!value?.value) return;
  const model = parseHelmDuplex(String(value.value));

  // Probe layout at the reference size, then crop the canvas to the natural
  // duplex extent so the resulting image isn't padded with blank pixels.
  const refLayout = drawDuplex(
    null as unknown as CanvasRenderingContext2D, 0, 0, REF_W, REF_H, model, undefined, true);
  const width = Math.max(REF_H, computeMaxChipEndX(refLayout) + RIGHT_PAD);

  const canvas = document.createElement('canvas');
  canvas.width = Math.max(1, Math.round((width + RIGHT_PAD) * IMAGE_SCALE));
  canvas.height = Math.max(1, Math.round(REF_H * IMAGE_SCALE));
  const g = canvas.getContext('2d');
  if (!g) return;
  g.scale(IMAGE_SCALE, IMAGE_SCALE);
  drawDuplex(g, 0, 0, width, REF_H, model);

  canvas.toBlob((blob) => {
    if (!blob) return;
    navigator.clipboard.write([new ClipboardItem({'image/png': blob})])
      .then(() => grok.shell.info('Image copied to clipboard'))
      .catch((er) => grok.shell.error(`Failed to copy image: ${er?.message ?? er}`));
  }, 'image/png');
}

/* -------------------------------------------------------------------------- *
 * Internals
 * -------------------------------------------------------------------------- */

function computeMaxChipEndX(layout: DuplexLayout): number {
  const lastSense = layout.senseChips.length ?
    layout.senseChips[layout.senseChips.length - 1] : null;
  const lastAnti = layout.antiChips.length ?
    layout.antiChips[layout.antiChips.length - 1] : null;
  return Math.max(
    (lastSense?.x ?? 0) + (lastSense?.w ?? 0),
    (lastAnti?.x ?? 0) + (lastAnti?.w ?? 0),
  );
}
