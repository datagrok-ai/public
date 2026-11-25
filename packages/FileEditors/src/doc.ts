import * as CFB from 'cfb';

/** Extracts the paragraphs of the main text of a Word 97-2003 binary document (.doc) from its piece table.
 * Formatting and images are dropped; a table row becomes one tab-separated paragraph.
 * See [MS-DOC] 2.4.1 "Retrieving Text". */
export function readDocParagraphs(bytes: Uint8Array): string[] {
  const cfb = CFB.read(bytes, {type: 'array'});
  const stream = (name: string) => {
    const content = CFB.find(cfb, `/${name}`)?.content;
    return content ? new Uint8Array(content as ArrayLike<number>) : null;
  };
  const word = stream('WordDocument');
  if (!word)
    throw new Error('Not a Word 97-2003 document: the WordDocument stream is missing');
  const fib = new DataView(word.buffer);
  if (fib.getUint16(0, true) !== 0xA5EC)
    throw new Error('Not a Word 97-2003 document');
  if (fib.getUint16(2, true) < 0xC1)
    throw new Error('Word 6/95 documents are not supported');
  const flags = fib.getUint16(0x0A, true);
  if (flags & 0x0100)
    throw new Error('Encrypted documents are not supported');
  const table = stream(flags & 0x0200 ? '1Table' : '0Table');
  if (!table)
    throw new Error('The table stream is missing');
  const ccpText = fib.getInt32(0x4C, true);
  const clx = new DataView(table.buffer, fib.getUint32(0x01A2, true), fib.getUint32(0x01A6, true));

  let pos = 0;
  while (clx.getUint8(pos) === 0x01)
    pos += 3 + clx.getUint16(pos + 1, true);
  if (clx.getUint8(pos) !== 0x02)
    throw new Error('The piece table is corrupted');
  const plc = pos + 5;
  const pieces = (clx.getUint32(pos + 1, true) - 4) / 12;
  const ansi = new TextDecoder('windows-1252');
  const utf16 = new TextDecoder('utf-16le');
  let text = '';
  for (let i = 0; i < pieces && text.length < ccpText; i++) {
    const cpStart = clx.getUint32(plc + i * 4, true);
    const count = Math.min(clx.getUint32(plc + (i + 1) * 4, true) - cpStart, ccpText - text.length);
    const fc = clx.getUint32(plc + (pieces + 1) * 4 + i * 8 + 2, true);
    const compressed = (fc & 0x40000000) !== 0;
    const offset = compressed ? (fc & 0x3FFFFFFF) / 2 : fc;
    text += compressed ? ansi.decode(word.subarray(offset, offset + count)) :
      utf16.decode(word.subarray(offset, offset + count * 2));
  }

  // A field is 0x13 code [0x14 result] 0x15 and may nest; only the results are text.
  const inCode: boolean[] = [];
  let visible = '';
  for (const c of text) {
    if (c === '\x13')
      inCode.push(true);
    else if (c === '\x14' && inCode.length)
      inCode[inCode.length - 1] = false;
    else if (c === '\x15')
      inCode.pop();
    else if (!inCode.includes(true))
      visible += c;
  }
  return visible
    .replace(/[\x00-\x06\x08\x1f]/g, '')
    .replace(/\x07\x07/g, '\r')
    .replace(/\x07/g, '\t')
    .replace(/[\x0b\x0c\x0e]/g, '\r')
    .replace(/\x1e/g, '-')
    .replace(/\r$/, '')
    .split('\r');
}
