import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {_package} from '../package';
import {parseFcs, fcsToDataFrame} from '../fcs/fcs-parser';


async function loadFcsBuffer(fileName: string): Promise<ArrayBuffer> {
  const bytes = await _package.files.readAsBytes(fileName);
  return bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer;
}


category('FCS Parser: Header', () => {
  test('parse header from 001.fcs', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    expect(parsed.header.version.startsWith('FCS'), true);
    expect(parsed.header.textStart > 0, true);
    expect(parsed.header.textEnd > parsed.header.textStart, true);
  });
});


category('FCS Parser: TEXT', () => {
  test('parse TEXT keywords from 001.fcs', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    expect(parsed.eventCount > 0, true);
    expect(parsed.paramCount > 0, true);
    expect(parsed.dataType, 'F');
  });

  test('parameter names extracted', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    expect(parsed.parameters.length, parsed.paramCount);
    for (const p of parsed.parameters)
      expect(p.shortName.length > 0, true);
  });
});


category('FCS Parser: DATA', () => {
  test('parse float data from 001.fcs', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    expect(parsed.data.length, parsed.paramCount);
    for (const col of parsed.data) {
      expect(col.length, parsed.eventCount);
      // Verify values are finite (not NaN or Infinity from bad parsing)
      let finiteCount = 0;
      for (let i = 0; i < Math.min(100, col.length); i++) {
        if (isFinite(col[i]))
          finiteCount++;
      }
      expect(finiteCount > 0, true);
    }
  });
});


category('FCS Import', () => {
  test('fcsToDataFrame creates correct DataFrame', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    const df = fcsToDataFrame(parsed);
    expect(df.rowCount, parsed.eventCount);
    expect(df.columns.length, parsed.paramCount);
    for (const p of parsed.parameters)
      expect(df.col(p.shortName) !== null, true);
  });

  test('DataFrame tags contain FCS keywords', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    const df = fcsToDataFrame(parsed);
    expect(df.getTag('$TOT'), String(parsed.eventCount));
    expect(df.getTag('$PAR'), String(parsed.paramCount));
  });

  test('column marker tags set from $PnS', async () => {
    const buffer = await loadFcsBuffer('001.fcs');
    const parsed = parseFcs(buffer);
    const df = fcsToDataFrame(parsed);
    for (const p of parsed.parameters) {
      const col = df.col(p.shortName)!;
      if (p.longName)
        expect(col.getTag('marker'), p.longName);
    }
  });
});
