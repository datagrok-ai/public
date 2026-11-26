import {roleOf} from './study-reader';

/** Names written differently in a group or animal column that are read as one (`as`). PowerPack turns a column of
 *  numbers into numbers, so "01" and "1" arrive as one: the workbook (a zip of XML files) is read again here as
 *  written, with what the browser has (DecompressionStream, DOMParser); anything it cannot read gives no finding. */
export interface Merged {role: 'group' | 'animal', texts: string[], as: string}

const u16 = (b: Uint8Array, i: number) => b[i] | b[i + 1] << 8;
const u32 = (b: Uint8Array, i: number) => (u16(b, i) | u16(b, i + 2) << 16) >>> 0;

/** The files of the zip that `wanted` names, as text (stored or deflated; anything else throws). */
export async function unzip(bytes: Uint8Array, wanted: (name: string) => boolean): Promise<Map<string, string>> {
  let end = bytes.length - 22;
  while (end >= 0 && u32(bytes, end) !== 0x06054b50)
    end--;
  if (end < 0)
    throw new Error('not a zip');
  const dec = new TextDecoder();
  const out = new Map<string, string>();
  let p = u32(bytes, end + 16);
  for (let k = u16(bytes, end + 10); k > 0; k--) {
    if (u32(bytes, p) !== 0x02014b50)
      throw new Error('the zip index does not read');
    const [method, size, nameLength, local] = [u16(bytes, p + 10), u32(bytes, p + 20), u16(bytes, p + 28),
      u32(bytes, p + 42)];
    const name = dec.decode(bytes.subarray(p + 46, p + 46 + nameLength));
    p += 46 + nameLength + u16(bytes, p + 30) + u16(bytes, p + 32);
    if (!wanted(name))
      continue;
    const start = local + 30 + u16(bytes, local + 26) + u16(bytes, local + 28);
    const data = bytes.subarray(start, start + size);
    if (method === 0)
      out.set(name, dec.decode(data));
    else if (method === 8) {
      out.set(name, await new Response(new Blob([data.slice()]).stream()
        .pipeThrough(new DecompressionStream('deflate-raw'))).text());
    } else
      throw new Error(`zip method ${method}`);
  }
  return out;
}

const xml = (text: string | undefined) => new DOMParser().parseFromString(text ?? '<x/>', 'application/xml');
/** The elements named `tag`, whatever prefix the writer gave the namespace. */
const all = (node: Document | Element, tag: string) => Array.from(node.getElementsByTagNameNS('*', tag));
/** A rich or plain text: its runs, without the phonetic guides some Japanese sheets add. */
const textOf = (el: Element) => all(el, 't').filter((t) => t.parentElement?.localName !== 'rPh')
  .map((t) => t.textContent ?? '').join('');
const columnOf = (ref: string) => [...ref.replace(/\d+$/, '')].reduce((n, c) => n * 26 + c.charCodeAt(0) - 64, 0) - 1;

/** The sheet's cells as written (the text of a text cell, a number as Excel stores it), row by row. */
export async function rawSheet(bytes: Uint8Array, sheet: string): Promise<string[][]> {
  const meta = await unzip(bytes, (n) => n === 'xl/workbook.xml' || n === 'xl/_rels/workbook.xml.rels' ||
    n === 'xl/sharedStrings.xml');
  // the sheet's link: its id attribute in the relationships namespace, whatever its prefix
  const id = Array.from(all(xml(meta.get('xl/workbook.xml')), 'sheet')
    .find((s) => s.getAttribute('name') === sheet)?.attributes ?? [])
    .find((a) => a.localName === 'id' && !!a.namespaceURI?.endsWith('/relationships'))?.value;
  const target = all(xml(meta.get('xl/_rels/workbook.xml.rels')), 'Relationship')
    .find((r) => r.getAttribute('Id') === id)?.getAttribute('Target');
  if (!target)
    throw new Error(`no sheet ${sheet}`);
  const path = target.startsWith('/') ? target.slice(1) : `xl/${target}`;
  const shared = all(xml(meta.get('xl/sharedStrings.xml')), 'si').map(textOf);
  return sheetCells((await unzip(bytes, (n) => n === path)).get(path) ?? '', shared);
}

/** A sheet's XML as its cells as written, row by row; `shared`, the workbook's shared texts. */
export function sheetCells(sheetXml: string, shared: string[]): string[][] {
  const rows: string[][] = [];
  // a row or a cell without its reference (some writers leave it out) is the one after the last
  let r = -1;
  for (const row of all(xml(sheetXml), 'row')) {
    r = row.hasAttribute('r') ? Number(row.getAttribute('r')) - 1 : r + 1;
    let col = -1;
    for (const c of all(row, 'c')) {
      col = c.hasAttribute('r') ? columnOf(c.getAttribute('r')!) : col + 1;
      if (!(r >= 0 && col >= 0))
        throw new Error(`cell ${c.getAttribute('r')}`);
      const t = c.getAttribute('t');
      const v = all(c, 'v')[0]?.textContent ?? '';
      (rows[r] ??= [])[col] = t === 's' ? shared[Number(v)] ?? '' : t === 'inlineStr' ? textOf(c) : v;
    }
  }
  return Array.from(rows, (row) => Array.from(row ?? [], (x) => (x ?? '').trim()));
}

// ponytail: a plain integer or decimal of up to 15 digits, a narrower set than any number Datagrok may read: a sign,
// "1e3" or ".5" give no finding, never a false one (and longer ids could round together in Number)
const NUMBER = /^-?\d{1,15}(\.\d{1,15})?$/;

/** In the group and animal columns (by their names, row 1), the names written differently that are read as one
 *  number: a column is read as numbers when every value is one. */
export function mergedNames(cells: string[][]): Merged[] {
  const out: Merged[] = [];
  (cells[0] ?? []).forEach((header, j) => {
    const role = roleOf(header)?.role;
    if (role !== 'group' && role !== 'animal')
      return;
    const values = cells.slice(1).map((row) => row[j] ?? '').filter((v) => v !== '');
    if (!values.length || !values.every((v) => NUMBER.test(v)))
      return;
    const by = new Map<string, Set<string>>();
    values.forEach((v) => by.set(String(Number(v)), (by.get(String(Number(v))) ?? new Set()).add(v)));
    by.forEach((texts, as) => texts.size > 1 && out.push({role, texts: [...texts], as}));
  });
  return out;
}

/** The merged names of the sheet `sheet` of the workbook; none when it cannot be read. */
export async function mergedInWorkbook(bytes: Uint8Array, sheet: string): Promise<Merged[]> {
  try {
    return mergedNames(await rawSheet(bytes, sheet));
  } catch (e) {
    console.warn('PMX Studio: the Excel sheet could not be read as written', e);
    return [];
  }
}

/** The warning line for the import summary: names, unless the groups are blinded; at most three. */
export function mergedWarning(merged: Merged[], blinded: boolean): string {
  if (!merged.length)
    return '';
  const quote = (m: Merged) => `${m.texts.map((t) => `"${t}"`).join(' and ')} (read as ${m.as})`;
  const shown = merged.filter((m) => !(blinded && m.role === 'group'));
  const hidden = merged.length - shown.length;
  const named = shown.slice(0, 3).map((m) => `${m.role} ${quote(m)}`);
  const more = shown.length - named.length + hidden;
  const list = named.length ?
    `${named.join('; ')}${more ? `; and ${more} more` : ''}` : `${hidden} group name${hidden > 1 ? 's' : ''}`;
  return `Names written differently in the sheet are read as one: ${list}. If it is a typo, correct it in the file; ` +
    'if they are different, give them names with letters.';
}
