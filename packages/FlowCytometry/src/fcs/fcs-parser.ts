import * as DG from 'datagrok-api/dg';

import {FcsHeader, FcsParameter, FcsSpillover, FcsParsed} from './fcs-types';


function readAscii(buffer: ArrayBuffer, start: number, length: number): string {
  return String.fromCharCode(...new Uint8Array(buffer, start, length));
}

function parseHeader(buffer: ArrayBuffer): FcsHeader {
  const version = readAscii(buffer, 0, 6).trim();
  if (!version.startsWith('FCS'))
    throw new Error(`Not an FCS file: expected "FCS" header, got "${version}"`);

  const field = (start: number) => {
    const s = readAscii(buffer, start, 8).trim();
    return s.length > 0 ? parseInt(s, 10) : 0;
  };

  return {
    version,
    textStart: field(10),
    textEnd: field(18),
    dataStart: field(26),
    dataEnd: field(34),
    analysisStart: field(42),
    analysisEnd: field(50),
  };
}

function parseTextSegment(buffer: ArrayBuffer, start: number, end: number): Map<string, string> {
  const bytes = new Uint8Array(buffer, start, end - start + 1);
  const text = String.fromCharCode(...bytes);
  const delimiter = text[0];
  const keywords = new Map<string, string>();

  // State-machine scanner to handle escaped delimiters (two consecutive = literal)
  const tokens: string[] = [];
  let current = '';
  let i = 1; // skip leading delimiter
  while (i < text.length) {
    if (text[i] === delimiter) {
      if (i + 1 < text.length && text[i + 1] === delimiter) {
        current += delimiter;
        i += 2;
      }
      else {
        tokens.push(current);
        current = '';
        i++;
      }
    }
    else {
      current += text[i];
      i++;
    }
  }
  if (current.length > 0)
    tokens.push(current);

  for (let t = 0; t + 1 < tokens.length; t += 2)
    keywords.set(tokens[t].toUpperCase(), tokens[t + 1]);

  return keywords;
}

function resolveDataOffsets(header: FcsHeader, keywords: Map<string, string>): {dataStart: number; dataEnd: number} {
  const beginData = keywords.get('$BEGINDATA');
  const endData = keywords.get('$ENDDATA');
  if (beginData && endData && (parseInt(beginData, 10) > 0 || header.dataStart === 0))
    return {dataStart: parseInt(beginData, 10), dataEnd: parseInt(endData, 10)};
  return {dataStart: header.dataStart, dataEnd: header.dataEnd};
}

function extractParameters(keywords: Map<string, string>, paramCount: number): FcsParameter[] {
  const params: FcsParameter[] = [];
  for (let i = 1; i <= paramCount; i++) {
    params.push({
      index: i,
      shortName: keywords.get(`$P${i}N`) ?? `P${i}`,
      longName: keywords.get(`$P${i}S`) ?? '',
      bits: parseInt(keywords.get(`$P${i}B`) ?? '32', 10),
      range: parseInt(keywords.get(`$P${i}R`) ?? '0', 10),
      amplification: keywords.get(`$P${i}E`) ?? '0,0',
    });
  }
  return params;
}

function parseDataSegment(
  buffer: ArrayBuffer, dataStart: number, dataEnd: number,
  eventCount: number, params: FcsParameter[],
  dataType: string, littleEndian: boolean
): Float32Array[] {
  const paramCount = params.length;
  const columns: Float32Array[] = [];
  for (let p = 0; p < paramCount; p++)
    columns.push(new Float32Array(eventCount));

  const view = new DataView(buffer, dataStart, dataEnd - dataStart + 1);
  let offset = 0;

  for (let event = 0; event < eventCount; event++) {
    for (let p = 0; p < paramCount; p++) {
      let value: number;
      if (dataType === 'F') {
        value = view.getFloat32(offset, littleEndian);
        offset += 4;
      }
      else if (dataType === 'D') {
        value = view.getFloat64(offset, littleEndian);
        offset += 8;
      }
      else {
        const bits = params[p].bits;
        if (bits <= 8) {
          value = view.getUint8(offset);
          offset += 1;
        }
        else if (bits <= 16) {
          value = view.getUint16(offset, littleEndian);
          offset += 2;
        }
        else if (bits <= 32) {
          value = view.getUint32(offset, littleEndian);
          offset += 4;
        }
        else {
          throw new Error(`Unsupported integer bit width: ${bits} for parameter ${params[p].shortName}`);
        }
        const range = params[p].range;
        if (range > 0 && (range & (range - 1)) !== 0)
          value = value & (range - 1);
      }
      columns[p][event] = value;
    }
  }

  return columns;
}

function parseSpillover(keywords: Map<string, string>): FcsSpillover | null {
  const raw = keywords.get('$SPILLOVER') ?? keywords.get('$SPILL') ?? keywords.get('SPILL');
  if (!raw)
    return null;

  const parts = raw.split(',').map((s) => s.trim());
  const n = parseInt(parts[0], 10);
  if (isNaN(n) || n <= 0)
    return null;

  const paramNames = parts.slice(1, 1 + n);
  const values = parts.slice(1 + n).map(Number);
  const matrix: number[][] = [];
  for (let row = 0; row < n; row++)
    matrix.push(values.slice(row * n, (row + 1) * n));

  return {n, paramNames, matrix};
}

export function parseFcs(buffer: ArrayBuffer): FcsParsed {
  const header = parseHeader(buffer);
  const keywords = parseTextSegment(buffer, header.textStart, header.textEnd);
  const {dataStart, dataEnd} = resolveDataOffsets(header, keywords);

  const eventCount = parseInt(keywords.get('$TOT') ?? '0', 10);
  const paramCount = parseInt(keywords.get('$PAR') ?? '0', 10);
  const dataType = (keywords.get('$DATATYPE') ?? 'F') as 'F' | 'D' | 'I';
  const byteOrd = keywords.get('$BYTEORD') ?? '1,2,3,4';
  const littleEndian = byteOrd.startsWith('1');

  const parameters = extractParameters(keywords, paramCount);
  const data = parseDataSegment(buffer, dataStart, dataEnd, eventCount, parameters, dataType, littleEndian);
  const spillover = parseSpillover(keywords);

  return {header, keywords, parameters, eventCount, paramCount, dataType, littleEndian, data, spillover};
}

export function fcsToDataFrame(parsed: FcsParsed): DG.DataFrame {
  const df = DG.DataFrame.create(parsed.eventCount);

  for (let i = 0; i < parsed.paramCount; i++) {
    const param = parsed.parameters[i];
    const col = DG.Column.fromFloat32Array(param.shortName, parsed.data[i]);
    if (param.longName)
      col.setTag('marker', param.longName);
    col.setTag('fcs.param.index', String(param.index));
    df.columns.add(col);
  }

  for (const [key, value] of parsed.keywords)
    df.setTag(key, value);

  if (parsed.spillover)
    df.setTag('spillover', JSON.stringify(parsed.spillover));

  df.name = parsed.keywords.get('$FIL') ?? parsed.keywords.get('FIL') ?? 'FCS Data';
  return df;
}
