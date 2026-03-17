import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';


const COLOR_SEM_TYPE = 'Color';

async function detect(values: (string | null)[]): Promise<string | null> {
  const col = DG.Column.fromStrings('color', values.map((v) => v ?? ''));
  return await grok.functions.call('PowerGrid:detectColorColumn', {col});
}

category('ColorCellRenderer', () => {
  test('detects hex6', async () => {
    const values = ['#FF5733', '#33FF57', '#3357FF', '#000000', '#FFFFFF', '#AABBCC', '#1A2B3C', '#FFA500', '#800080', '#DC143C'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects hex3', async () => {
    const values = ['#F00', '#0F0', '#00F', '#FFF', '#000', '#ABC', '#369', '#C0C', '#F80', '#1E1'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects hex8 with alpha', async () => {
    const values = ['#FF000080', '#00FF00FF', '#0000FF40', '#FFFFFF00', '#00000099', '#AABBCCDD', '#11223344', '#DEADBEEF', '#CAFEBABE', '#12345678'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects rgb()', async () => {
    const values = ['rgb(255, 0, 0)', 'rgb(0, 255, 0)', 'rgb(0, 0, 255)', 'rgb(128, 128, 128)', 'rgb(255, 165, 0)',
      'rgb(75, 0, 130)', 'rgb(0, 139, 139)', 'rgb(210, 180, 140)', 'rgb(100%, 0%, 50%)', 'rgb(64, 224, 208)'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects rgba()', async () => {
    const values = ['rgba(255, 0, 0, 0.5)', 'rgba(0, 128, 255, 0.8)', 'rgba(0, 0, 0, 0)', 'rgba(255, 255, 255, 1)',
      'rgba(34, 139, 34, 0.75)', 'rgba(147, 112, 219, 0.85)', 'rgba(255, 215, 0, 1)',
      'rgba(0, 0, 0, 0.3)', 'rgba(100, 100, 100, 0.5)', 'rgba(200, 50, 50, 0.9)'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects hsl()', async () => {
    const values = ['hsl(0, 100%, 50%)', 'hsl(120, 100%, 50%)', 'hsl(240, 100%, 50%)', 'hsl(60, 100%, 50%)', 'hsl(300, 100%, 25%)',
      'hsl(180, 50%, 70%)', 'hsl(30, 80%, 60%)', 'hsl(210, 100%, 56%)', 'hsl(275, 60%, 45%)', 'hsl(15, 80%, 55%)'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects hsla()', async () => {
    const values = ['hsla(0, 100%, 50%, 0.5)', 'hsla(240, 100%, 50%, 0.3)', 'hsla(120, 60%, 40%, 0.9)', 'hsla(180, 50%, 50%, 0.7)',
      'hsla(30, 80%, 60%, 1)', 'hsla(300, 100%, 25%, 0.4)', 'hsla(60, 100%, 50%, 0.8)', 'hsla(210, 100%, 56%, 0.6)',
      'hsla(275, 60%, 45%, 0.2)', 'hsla(15, 80%, 55%, 0.5)'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects mixed formats', async () => {
    const values = ['#FF5733', 'rgb(0, 255, 0)', 'hsl(240, 100%, 50%)', '#FFF', 'rgba(0, 0, 0, 0.5)',
      '#AABBCC', 'hsla(120, 60%, 40%, 0.9)', '#00FF00FF', 'rgb(128, 128, 128)', 'hsl(60, 100%, 50%)'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('detects with some empty values', async () => {
    const values = ['#FF5733', '', '#33FF57', '', '#3357FF', '#000000', '', '#FFFFFF', '#AABBCC', '#1A2B3C', '#FFA500',
      '#800080', '#DC143C', '#7FFFD4', '#FF69B4'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('accepts out-of-range values that browsers clamp', async () => {
    const values = ['rgb(256, 0, 0)', 'rgb(-10, 200, 50)', 'hsl(720, 100%, 50%)', 'rgba(0, 0, 0, 2)',
      '#FF5733', '#33FF57', '#3357FF', '#000000', '#FFFFFF', '#AABBCC'];
    expect(await detect(values), COLOR_SEM_TYPE);
  });

  test('does not detect with invalid values mixed in', async () => {
    // ratio=1 by default — any sampled invalid value stops detection
    const values = ['#FF5733', 'notacolor', '#33FF57', '#3357FF', 'COLOR', '#000000', '#FFFFFF', '#AABBCC', '255 0 0',
      '#1A2B3C', '#FFA500', '#800080', '#DC143C', '#7FFFD4', '#FF69B4', '#4B0082', '#FF4500', '#2E8B57', '#F0F0F0', '#123456'];
    expect(await detect(values), null);
  });

  test('does not detect named CSS colors', async () => {
    const values = ['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'teal', 'pink', 'lime', 'navy'];
    expect(await detect(values), null);
  });

  test('does not detect plain text', async () => {
    const values = ['apple', 'banana', 'cherry', 'date', 'elderberry', 'fig', 'grape', 'honeydew', 'kiwi', 'lemon'];
    expect(await detect(values), null);
  });

  test('does not detect numeric strings', async () => {
    const values = ['123', '456', '789', '1011', '1213', '1415', '1617', '1819', '2021', '2223'];
    expect(await detect(values), null);
  });

  test('does not detect invalid hex strings', async () => {
    const values = ['#FE', '#GGGGGG', '#12345', '#1234567', '#GGG', '#XYZ', '#ZZZZZZ', '#!!@@##', '#12G', '#AABBCX'];
    expect(await detect(values), null);
  });

  test('does not detect with only one distinct value', async () => {
    const values = ['notacolor', 'notacolor', 'notacolor', 'notacolor', 'notacolor'];
    expect(await detect(values), null);
  });

  test('detects color column in color-test.csv', async () => {
    const csv = await _package.files.readAsText('color-test.csv');
    const df = DG.DataFrame.fromCsv(csv);
    await grok.data.detectSemanticTypes(df);
    expect(df.getCol('color').semType, COLOR_SEM_TYPE);
  });
});
