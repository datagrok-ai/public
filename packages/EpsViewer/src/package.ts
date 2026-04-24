import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import UTIF from 'utif2';

import '../css/eps-viewer.css';

export * from './package.g';
export const _package = new DG.Package();

declare const FromPS: any;
declare const ToContext2D: any;
declare const ToPDF: any;

const DOS_EPS_MAGIC = [0xc5, 0xd0, 0xd3, 0xc6];

let udocLoaded = false;

async function loadUdoc(): Promise<void> {
  if (udocLoaded)
    return;
  await DG.Utils.loadJsCss([
    `${_package.webRoot}/dist/UDOC.js`,
    `${_package.webRoot}/dist/FromPS.js`,
    `${_package.webRoot}/dist/ToContext2D.js`,
  ]);
  if (typeof ToPDF === 'undefined')
    (window as any).ToPDF = {_flt: (n: number) => '' + parseFloat(n.toFixed(2))};
  udocLoaded = true;
}

function readU32LE(bytes: Uint8Array, off: number): number {
  return bytes[off] | (bytes[off + 1] << 8) | (bytes[off + 2] << 16) | (bytes[off + 3] << 24);
}

/// DOS EPS binary files start with C5 D0 D3 C6 and embed a TIFF or WMF preview.
/// Header: 4-byte magic, then uint32 LE: psOffset, psLen, wmfOffset, wmfLen, tiffOffset, tiffLen.
function extractDosEpsTiff(bytes: Uint8Array): Uint8Array | null {
  if (bytes.length < 32)
    return null;
  for (let i = 0; i < 4; i++)
    if (bytes[i] !== DOS_EPS_MAGIC[i])
      return null;
  const tiffOffset = readU32LE(bytes, 20);
  const tiffLen = readU32LE(bytes, 24);
  if (tiffOffset === 0 || tiffLen === 0 || tiffOffset + tiffLen > bytes.length)
    return null;
  return bytes.subarray(tiffOffset, tiffOffset + tiffLen);
}

/// ASCII EPS files can embed a hex-encoded preview between %%BeginPreview and %%EndPreview.
function extractEpsiPreview(bytes: Uint8Array): HTMLCanvasElement | null {
  const header = new TextDecoder('latin1').decode(bytes.subarray(0, Math.min(bytes.length, 8192)));
  const begin = header.indexOf('%%BeginPreview:');
  if (begin < 0)
    return null;
  const lineEnd = header.indexOf('\n', begin);
  const params = header.substring(begin + 15, lineEnd).trim().split(/\s+/).map(Number);
  if (params.length < 4 || params.some(isNaN))
    return null;
  const [w, h, depth] = params;
  const fullText = new TextDecoder('latin1').decode(bytes);
  const dataStart = fullText.indexOf('\n', begin) + 1;
  const dataEnd = fullText.indexOf('%%EndPreview', dataStart);
  if (dataEnd < 0)
    return null;
  const hex = fullText.substring(dataStart, dataEnd).replace(/[^0-9a-fA-F]/g, '');
  const bitsPerPixel = depth;
  const rowBytes = Math.ceil(w * bitsPerPixel / 8);
  const pixelData = new Uint8Array(rowBytes * h);
  for (let i = 0; i < pixelData.length && i * 2 + 1 < hex.length; i++)
    pixelData[i] = parseInt(hex.substring(i * 2, i * 2 + 2), 16);
  const canvas = document.createElement('canvas');
  canvas.width = w;
  canvas.height = h;
  const ctx = canvas.getContext('2d')!;
  const img = ctx.createImageData(w, h);
  for (let y = 0; y < h; y++) {
    for (let x = 0; x < w; x++) {
      let v = 0;
      if (bitsPerPixel === 1) {
        const byte = pixelData[y * rowBytes + (x >> 3)];
        v = ((byte >> (7 - (x & 7))) & 1) ? 255 : 0;
      } else if (bitsPerPixel === 8)
        v = pixelData[y * rowBytes + x];
      const p = (y * w + x) * 4;
      img.data[p] = img.data[p + 1] = img.data[p + 2] = v;
      img.data[p + 3] = 255;
    }
  }
  ctx.putImageData(img, 0, 0);
  return canvas;
}

function renderTiff(tiff: Uint8Array): HTMLCanvasElement {
  const slice = new Uint8Array(tiff).buffer;
  const ifds = UTIF.decode(slice);
  UTIF.decodeImage(slice, ifds[0]);
  const rgba = UTIF.toRGBA8(ifds[0]);
  const canvas = document.createElement('canvas');
  canvas.width = ifds[0].width;
  canvas.height = ifds[0].height;
  const ctx = canvas.getContext('2d')!;
  const img = new ImageData(new Uint8ClampedArray(rgba), ifds[0].width, ifds[0].height);
  ctx.putImageData(img, 0, 0);
  return canvas;
}

async function renderEps(bytes: Uint8Array): Promise<HTMLCanvasElement> {
  const dosTiff = extractDosEpsTiff(bytes);
  if (dosTiff)
    return renderTiff(dosTiff);
  const epsiPreview = extractEpsiPreview(bytes);
  if (epsiPreview)
    return epsiPreview;
  await loadUdoc();
  const renderer = new ToContext2D(0, 1);
  FromPS.Parse(bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength), renderer);
  return renderer.canvas as HTMLCanvasElement;
}

export class PackageFunctions {
  @grok.decorators.fileViewer({fileViewer: 'eps'})
  static previewEps(file: DG.FileInfo): DG.View {
    const view = DG.View.create();
    view.name = file.name;

    const host = ui.div([], 'eps-viewer-host');
    view.append(host);

    file.readAsBytes().then(async (bytes) => {
      try {
        const canvas = await renderEps(bytes);
        canvas.removeAttribute('style');
        host.appendChild(canvas);
      }
      catch (err) {
        const msg = err instanceof Error ? err.message : String(err);
        grok.shell.error(`EPS preview error: ${msg}`);
        host.appendChild(ui.divText('Failed to render EPS file.', 'eps-viewer-error'));
      }
    });

    return view;
  }
}
