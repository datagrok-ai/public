/* A Word document (.docx) from a few blocks: headings, paragraphs, tables and checkboxes, with the fields left to the
   user as grey "[…]" runs that Word's search finds. Written by hand, as xlsx-raw.ts reads: the XML parts in a zip,
   with each part's CRC-32. No platform imports. */

/** A run of text: plain, bold, or a field the user fills (grey, in brackets). A '\n' breaks the line. */
import {maxOf} from './format';

export type Run = string | {text: string, bold?: boolean, fill?: boolean};
/** A table cell: one run or several. */
export type Cell = Run | Run[];
export type Block =
  | {heading: string, level?: 1 | 2 | 3}
  | {p: Cell, small?: boolean}
  | {check: boolean, text: Cell}
  /** the first row is the header when `head` (repeated on each page); `widths` relative, one per column */
  | {rows: Cell[][], head?: boolean, widths?: number[]};

/** A field the user fills: '[Grant code]', grey. */
export const fill = (label: string): Run => ({text: `[${label}]`, fill: true});

const esc = (s: string) => s.replace(/[\u0000-\u0008\u000b\u000c\u000e-\u001f]/g, '').replace(/&/g, '&amp;')
  .replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

const list = (c: Cell): Run[] => Array.isArray(c) ? c : [c];
const runs = (c: Cell): string => list(c).map((r) => {
  const x = typeof r === 'string' ? {text: r} : r;
  const props = (x.bold ? '<w:b/>' : '') + (x.fill ? '<w:color w:val="808080"/>' : '');
  return x.text.split('\n').map((t, i) => `<w:r>${props ? `<w:rPr>${props}</w:rPr>` : ''}${i ? '<w:br/>' : ''}` +
    `<w:t xml:space="preserve">${esc(t)}</w:t></w:r>`).join('');
}).join('');

const para = (c: Cell, style?: string) => `<w:p>${style ? `<w:pPr><w:pStyle w:val="${style}"/></w:pPr>` : ''}` +
  `${runs(c)}</w:p>`;
const bold = (c: Cell): Run[] => list(c).map((x) => typeof x === 'string' ? {text: x, bold: true} : {...x, bold: true});

function table(b: {rows: Cell[][], head?: boolean, widths?: number[]}, width: number): string {
  const cols = maxOf(b.rows.map((r) => r.length));
  const rel = b.widths?.length === cols ? b.widths : Array(cols).fill(1);
  const sum = rel.reduce((a, x) => a + x, 0);
  const tw = rel.map((x) => Math.floor(width * x / sum));
  const row = (r: Cell[], i: number) => {
    const head = !!b.head && i === 0;
    // a row stays on one page
    return `<w:tr><w:trPr><w:cantSplit/>${head ? '<w:tblHeader/>' : ''}</w:trPr>${tw.map((w, j) =>
      `<w:tc><w:tcPr><w:tcW w:w="${w}" w:type="dxa"/>${head ? '<w:shd w:val="clear" w:color="auto" ' +
        'w:fill="F0F0F0"/>' : ''}</w:tcPr>${para(head ? bold(r[j] ?? '') : r[j] ?? '')}</w:tc>`).join('')}</w:tr>`;
  };
  // a table needs a paragraph after it, or Word joins it with the next one
  return `<w:tbl><w:tblPr><w:tblStyle w:val="TableGrid"/><w:tblW w:w="${width}" w:type="dxa"/>` +
    `<w:tblLayout w:type="fixed"/></w:tblPr><w:tblGrid>${tw.map((w) => `<w:gridCol w:w="${w}"/>`).join('')}` +
    `</w:tblGrid>${b.rows.map(row).join('')}</w:tbl><w:p/>`;
}

const NS = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"';
// A4 in twentieths of a point, 2 cm margins
const [A4W, A4H, MARGIN] = [11906, 16838, 1134];

/** The document's main part: the blocks, then the page (A4, `landscape` for wide tables). */
export function documentXml(blocks: Block[], landscape = false): string {
  const [w, h] = landscape ? [A4H, A4W] : [A4W, A4H];
  const body = blocks.map((b) => 'heading' in b ? para(b.heading, `Heading${b.level ?? 1}`) :
    'p' in b ? para(b.p, b.small ? 'Small' : undefined) :
      'check' in b ? para([b.check ? '☒ ' : '☐ ', ...list(b.text)]) : table(b, w - 2 * MARGIN)).join('');
  return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document ${NS}><w:body>${body}` +
    `<w:sectPr><w:pgSz w:w="${w}" w:h="${h}"${landscape ? ' w:orient="landscape"' : ''}/>` +
    `<w:pgMar w:top="${MARGIN}" w:right="${MARGIN}" w:bottom="${MARGIN}" w:left="${MARGIN}" w:header="709" ` +
    'w:footer="709" w:gutter="0"/></w:sectPr></w:body></w:document>';
}

const style = (id: string, name: string, rpr: string, ppr = '') => `<w:style w:type="paragraph" w:styleId="${id}">` +
  `<w:name w:val="${name}"/><w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:qFormat/>${ppr}` +
  `<w:rPr>${rpr}</w:rPr></w:style>`;
const heading = (n: number, size: number, before: number) => style(`Heading${n}`, `heading ${n}`,
  `<w:b/><w:sz w:val="${size}"/>`, `<w:pPr><w:keepNext/><w:spacing w:before="${before}" w:after="80"/>` +
  `<w:outlineLvl w:val="${n - 1}"/></w:pPr>`);
const STYLES = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:styles ${NS}><w:docDefaults><w:rPrDefault>` +
  '<w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri" w:eastAsia="Calibri"/><w:sz w:val="20"/>' +
  '<w:lang w:val="en-GB"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="80"/></w:pPr>' +
  '</w:pPrDefault></w:docDefaults>' +
  '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>' +
  heading(1, 30, 240) + heading(2, 24, 200) + heading(3, 20, 160) +
  style('Small', 'Small', '<w:color w:val="666666"/><w:sz w:val="16"/>') +
  '<w:style w:type="table" w:styleId="TableGrid"><w:name w:val="Table Grid"/><w:tblPr><w:tblBorders>' +
  ['top', 'left', 'bottom', 'right', 'insideH', 'insideV'].map((s) =>
    `<w:${s} w:val="single" w:sz="4" w:space="0" w:color="A0A0A0"/>`).join('') +
  '</w:tblBorders><w:tblCellMar><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tblCellMar>' +
  '</w:tblPr></w:style></w:styles>';

const REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
const PKG = 'http://schemas.openxmlformats.org/package/2006';
const WML = 'application/vnd.openxmlformats-officedocument.wordprocessingml';
const XML = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
export const DOCX_TYPE = `${WML}.document`;

/** The .docx file of the blocks, `title` in the document's properties. */
export function docx(blocks: Block[], o: {title?: string, landscape?: boolean} = {}): Promise<Uint8Array> {
  return zip([
    ['[Content_Types].xml', `${XML}<Types xmlns="${PKG}/content-types">` +
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' +
      '<Default Extension="xml" ContentType="application/xml"/>' +
      `<Override PartName="/word/document.xml" ContentType="${WML}.document.main+xml"/>` +
      `<Override PartName="/word/styles.xml" ContentType="${WML}.styles+xml"/>` +
      '<Override PartName="/docProps/core.xml" ' +
      'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/></Types>'],
    ['_rels/.rels', `${XML}<Relationships xmlns="${PKG}/relationships">` +
      `<Relationship Id="rId1" Type="${REL}/officeDocument" Target="word/document.xml"/>` +
      `<Relationship Id="rId2" Type="${PKG}/relationships/metadata/core-properties" Target="docProps/core.xml"/>` +
      '</Relationships>'],
    ['docProps/core.xml', `${XML}<cp:coreProperties xmlns:cp="${PKG}/metadata/core-properties" ` +
      `xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>${esc(o.title ?? '')}</dc:title>` +
      '<dc:creator>PMX Studio</dc:creator></cp:coreProperties>'],
    ['word/document.xml', documentXml(blocks, o.landscape)],
    ['word/styles.xml', STYLES],
    ['word/_rels/document.xml.rels', `${XML}<Relationships xmlns="${PKG}/relationships">` +
      `<Relationship Id="rId1" Type="${REL}/styles" Target="styles.xml"/></Relationships>`],
  ]);
}

const CRC = Array.from({length: 256}, (_, n) => {
  for (let k = 0; k < 8; k++)
    n = n & 1 ? 0xedb88320 ^ (n >>> 1) : n >>> 1;
  return n >>> 0;
});
export function crc32(b: Uint8Array): number {
  let c = 0xffffffff;
  for (const x of b)
    c = CRC[(c ^ x) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

/** Little-endian fields [offset, value, bytes] into a header of `size` bytes, `name` after them. */
function header(size: number, fields: number[][], name: Uint8Array): Uint8Array {
  const h = new Uint8Array(size + name.length);
  const v = new DataView(h.buffer);
  fields.forEach(([at, x, s]) => s === 4 ? v.setUint32(at, x >>> 0, true) : v.setUint16(at, x, true));
  h.set(name, size);
  return h;
}

const deflate = async (b: Uint8Array) => new Uint8Array(await new Response(new Blob([b.slice()]).stream()
  .pipeThrough(new CompressionStream('deflate-raw'))).arrayBuffer());

/** A zip of the files (name, text), names in UTF-8, each deflated with what the browser has (as xlsx-raw.ts
 *  inflates); dated 1980-01-01, the earliest DOS date, so the same document gives the same bytes. */
export async function zip(files: [string, string][]): Promise<Uint8Array> {
  const enc = new TextEncoder();
  const parts: Uint8Array[] = [];
  const central: Uint8Array[] = [];
  let offset = 0;
  const DATE = 0x21;
  for (const [n, text] of files) {
    const name = enc.encode(n);
    const raw = enc.encode(text);
    const data = await deflate(raw);
    const crc = crc32(raw);
    // version 2.0, flag 0x800 (UTF-8 names), method 8 (deflate), time 0, compressed then full size
    const local = header(30, [[0, 0x04034b50, 4], [4, 20, 2], [6, 0x800, 2], [8, 8, 2], [12, DATE, 2],
      [14, crc, 4], [18, data.length, 4], [22, raw.length, 4], [26, name.length, 2]], name);
    central.push(header(46, [[0, 0x02014b50, 4], [4, 20, 2], [6, 20, 2], [8, 0x800, 2], [10, 8, 2],
      [14, DATE, 2], [16, crc, 4], [20, data.length, 4], [24, raw.length, 4], [28, name.length, 2],
      [42, offset, 4]], name));
    parts.push(local, data);
    offset += local.length + data.length;
  }
  const size = central.reduce((a, x) => a + x.length, 0);
  const end = header(22, [[0, 0x06054b50, 4], [8, files.length, 2], [10, files.length, 2], [12, size, 4],
    [16, offset, 4]], new Uint8Array());
  const out = new Uint8Array(offset + size + end.length);
  let p = 0;
  for (const x of [...parts, ...central, end]) {
    out.set(x, p);
    p += x.length;
  }
  return out;
}
