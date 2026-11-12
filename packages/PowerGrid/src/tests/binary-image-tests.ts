import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {waitUntil, category, expect, expectArray, test} from '@datagrok-libraries/test/src/test';
import {testData} from '@datagrok-libraries/test/src/dataframe-utils';
import {BinaryImageCellRenderer} from '../cell-types/binary-image-cell-renderer';
import {projectRoundTrip} from './utils';

const ICO =
  'AAABAAEAEAsAAAAAIAAtAQAAFgAAAIlQTkcNChoKAAAADUlIRFIAAAAQAAAACwgCAAAA+YCabgAAAPRJREFUeJxjlFuwioGB4T8D' +
  'YcAIJlkgHCYGhr/////7j10jEyMjMyMjRI7l7///PKwse/09xLk4WJiYsGr48+/fy28/nDfu+PL7D8iGf///S/Nwbbj36MyrNzys' +
  'LH+RrGFmZPjy+4+JmEiAkhzEfogGkNzauw+XXLvFx8X+B0kHCzPjp28/Y7TUwBqQ/MDAwCDCyS7OyyPCyfbnH5IGJsY3zKwinOwI' +
  'ETjrzfefLz9/+f4Xiw1vvv9E0cAEDrBgZXl5Xm5cfoArg2hgfPrlm4+CTICSHK5QevrlGxMjSAcjJOIYiY8HiNA/BgZGsChWDVAF' +
  'YAYAggpo7yDO8agAAAAASUVORK5CYII=';
const AVIF =
  'AAAAIGZ0eXBhdmlmAAAAAGF2aWZtaWYxbWlhZk1BMUIAAADrbWV0YQAAAAAAAAAhaGRscgAAAAAAAAAAcGljdAAAAAAAAAAAAAAA' +
  'AAAAAAAOcGl0bQAAAAAAAQAAAB5pbG9jAAAAAEQAAAEAAQAAAAEAAAETAAAANwAAAChpaW5mAAAAAAABAAAAGmluZmUCAAAAAAEA' +
  'AGF2MDFDb2xvcgAAAABqaXBycAAAAEtpcGNvAAAAFGlzcGUAAAAAAAAAGAAAABAAAAAQcGl4aQAAAAADCAgIAAAADGF2MUOBAAwA' +
  'AAAAE2NvbHJuY2x4AAEADQAGgAAAABdpcG1hAAAAAAAAAAEAAQQBAoMEAAAAP21kYXQSAAoJGBDv7RAQ0GhAMihFMeHhmUCCCCeU' +
  'AAAAV98kVTsMVdTOGBM6Mnbv+0tkUC2Q7fBwTiPg';

const fromBase64 = (s: string) => Uint8Array.from(atob(s), (c) => c.charCodeAt(0));

async function redPng(): Promise<Uint8Array> {
  const canvas = new OffscreenCanvas(40, 20);
  const g = canvas.getContext('2d')!;
  g.fillStyle = 'rgb(255, 0, 0)';
  g.fillRect(0, 0, 40, 20);
  return new Uint8Array(await (await canvas.convertToBlob({type: 'image/png'})).arrayBuffer());
}

function bytesColumn(name: string, values: (Uint8Array | null)[]): DG.Column {
  return DG.Column.fromList(DG.COLUMN_TYPE.BYTE_ARRAY, name, values);
}

category('Binary image', () => {
  test('detects columns with images only', async () => {
    const png = await redPng();
    const pdf = new TextEncoder().encode('%PDF-1.4 1 0 obj endobj');
    const df = DG.DataFrame.fromColumns([
      bytesColumn('png', [png, null, null]),
      bytesColumn('ico', [fromBase64(ICO), null, null]),
      bytesColumn('avif', [fromBase64(AVIF), null, null]),
      bytesColumn('mixed', [pdf, null, png]),
      bytesColumn('empty', [null, null, null]),
      bytesColumn('pdf', [pdf, pdf, null]),
    ]);
    await grok.data.detectSemanticTypes(df);
    for (const name of ['png', 'ico', 'avif', 'mixed', 'empty'])
      expect(df.col(name)!.semType, 'BinaryImage', name);
    expect(df.col('pdf')!.semType === 'BinaryImage', false, 'pdf');
  });

  // the contract of the library's detector auto-tests, which skip detectors declared with meta.role
  test('tags only the BinaryImage column of the detectors test data', async () => {
    const tagged: string[] = [];
    for (const col of testData.clone().columns) {
      if (await grok.functions.call('PowerGrid:detectByteArrayColumn', {col}))
        tagged.push(col.name);
    }
    expectArray(tagged, ['BinaryImage']);
  });

  test('draws a decoded image in the same paint', async () => {
    const bytes = await redPng();
    const gridCell = DG.Grid.create(DG.DataFrame.fromColumns([bytesColumn('image', [bytes])])).cell('image', 0);
    await waitUntil(() => BinaryImageCellRenderer.getImage(bytes).naturalWidth > 0, 'image not decoded', 3000);
    const g = ui.canvas(100, 50).getContext('2d')!;
    new BinaryImageCellRenderer().render(g, 0, 0, 100, 50, gridCell, gridCell.style);
    expectArray(Array.from(g.getImageData(50, 25, 1, 1).data), [255, 0, 0, 255]);
  });

  test('draws into a canvas that is not repainted', async () => {
    const bytes = await redPng();
    const gridCell = DG.Grid.create(DG.DataFrame.fromColumns([bytesColumn('image', [bytes])])).cell('image', 0);
    const g = ui.canvas(100, 50).getContext('2d')!;
    new BinaryImageCellRenderer().render(g, 0, 0, 100, 50, gridCell, gridCell.style);
    expect(g.getImageData(50, 25, 1, 1).data[3], 0, 'drawn before the image was decoded');
    await waitUntil(() => g.getImageData(50, 25, 1, 1).data[0] === 255, 'the image was not drawn once loaded', 3000);
  });

  test('survives a project round trip', async () => {
    const png = await redPng();
    const view = grok.shell.addTableView(DG.DataFrame.fromColumns([
      DG.Column.fromStrings('name', ['red', 'none']), bytesColumn('image', [png, null])]));
    await projectRoundTrip(view, async (restored) => {
      const col = restored.dataFrame.col('image')!;
      // the project brings the table back before its column types
      await waitUntil(() => col.semType === 'BinaryImage', () => `the column came back as ${col.semType}`, 5000);
      expectArray(Array.from(col.get(0) as Uint8Array), Array.from(png));
      expect(col.isNone(1), true);
      await waitUntil(() => restored.grid.col('image')?.renderer?.cellType === 'BinaryImage',
        'the column is not drawn as images', 5000);
    });
  });

  test('keeps a bounded number of images', async () => {
    const max = BinaryImageCellRenderer.maxImages;
    const first = new Uint8Array(20);
    BinaryImageCellRenderer.getImage(first);
    for (let i = 0; i < max; i++)
      BinaryImageCellRenderer.getImage(new Uint8Array(20));
    expect(BinaryImageCellRenderer.images.size, max);
    expect(BinaryImageCellRenderer.images.has(first), false);
  });
});
