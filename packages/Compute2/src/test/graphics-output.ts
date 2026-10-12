import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type ExcelJS from 'exceljs';
import {category, test, expect, awaitCheck} from '@datagrok-libraries/test/src/test';
import {richFunctionViewReport} from '@datagrok-libraries/compute-utils';
import {closeViewsAfter, awaitWebComponents} from './utils';

const PNG_PREFIX = 'iVBORw0KGgo';

const plotValue = (call: DG.FuncCall) => String(call.outputs.get('plot') ?? '');

async function openEditor(call: DG.FuncCall): Promise<DG.ViewBase> {
  await awaitWebComponents();
  const view = await grok.functions.call('Compute2:RichFunctionViewEditor', {call}) as unknown as DG.ViewBase;
  await awaitCheck(() => view.root.querySelector('dg-input-form') !== null, 'RFV form not rendered', 15000);
  return view;
}

async function run(view: DG.ViewBase) {
  await awaitCheck(() => Array.from(view.root.querySelectorAll('button')).some((b) => b.textContent?.trim() === 'Run'),
    'Run button not found', 10000);
  Array.from(view.root.querySelectorAll('button')).find((b) => b.textContent?.trim() === 'Run')!.click();
}

const plotImg = (view: DG.ViewBase) =>
  view.root.querySelector('[dock-spawn-title="plot"] img') as HTMLImageElement | null;

category('RFV: graphics outputs', () => {
  const track = closeViewsAfter();

  test('PNG output renders as an image once the run fills it', async () => {
    const call = DG.Func.byName('Compute2:GraphicsOutputTest').prepare({a: 7});
    const view = track(await openEditor(call));
    expect(view.root.querySelector('[dock-spawn-title="plot"]') === null, true,
      'empty plot should be hidden before the run');

    await run(view);
    await awaitCheck(() => plotValue(call).startsWith(PNG_PREFIX), 'run produced no plot', 15000);
    await awaitCheck(() => plotImg(view!) !== null, 'plot image not rendered', 15000);
    expect(plotImg(view)!.src.startsWith(`data:image/png;base64,${PNG_PREFIX}`), true, 'PNG data URL expected');
    expect((view.root.textContent ?? '').includes(PNG_PREFIX), false, 'base64 text must not be shown');
  });

  test('SVG output renders as an SVG image', async () => {
    const call = DG.Func.byName('Compute2:GraphicsSvgOutputTest').prepare({a: 1});
    const view = track(await openEditor(call));
    await run(view);
    await awaitCheck(() => plotValue(call).startsWith('<svg'), 'run produced no plot', 15000);
    await awaitCheck(() => plotImg(view!) !== null, 'plot image not rendered', 15000);
    expect(plotImg(view)!.src.startsWith('data:image/svg+xml'), true, 'SVG data URL expected');
    expect((view.root.textContent ?? '').includes('<svg'), false, 'SVG markup must not be shown as text');
  });
});

// exceljs typings omit `ext`, which addImage({tl, ext}) keeps on the range
const imageSize = (range: unknown) => (range as {ext?: {width: number, height: number}}).ext;

async function exportWorkbook(call: DG.FuncCall): Promise<ExcelJS.Workbook> {
  const [, wb] = await richFunctionViewReport('Excel', call.func, call, {});
  return wb;
}

category('Export: graphics outputs', () => {
  test('PNG output gets its own sheet at native size', async () => {
    const call = await DG.Func.byName('Compute2:GraphicsOutputTest').prepare({a: 7}).call();
    expect(plotValue(call).startsWith(PNG_PREFIX), true, 'run produced no plot');
    const wb = await exportWorkbook(call);
    const sheet = wb.getWorksheet('plot');
    expect(sheet != null, true, 'plot sheet expected');
    const images = sheet!.getImages();
    expect(images.length, 1);
    expect(imageSize(images[0].range)?.width, 3);
    expect(imageSize(images[0].range)?.height, 2);
  });

  test('SVG output is fitted to the chart box keeping its aspect ratio', async () => {
    const call = await DG.Func.byName('Compute2:GraphicsSvgOutputTest').prepare({a: 1}).call();
    expect(plotValue(call).startsWith('<svg'), true, 'run produced no plot');
    const wb = await exportWorkbook(call);
    const images = wb.getWorksheet('plot')?.getImages() ?? [];
    expect(images.length, 1);
    expect(imageSize(images[0].range)?.width, 1280);
    expect(imageSize(images[0].range)?.height, 640);
  });

  test('Empty graphics output gets no sheet', async () => {
    const call = DG.Func.byName('Compute2:GraphicsOutputTest').prepare({a: 7});
    const wb = await exportWorkbook(call);
    expect(wb.getWorksheet('plot') == null, true, 'no plot sheet for an unrun call');
  });
});
