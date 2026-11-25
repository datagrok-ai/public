import JSZip from 'jszip';
import {isSafeUrl} from './docx';

const TEXT = 'urn:oasis:names:tc:opendocument:xmlns:text:1.0';
const TABLE = 'urn:oasis:names:tc:opendocument:xmlns:table:1.0';
const DRAW = 'urn:oasis:names:tc:opendocument:xmlns:drawing:1.0';
const STYLE = 'urn:oasis:names:tc:opendocument:xmlns:style:1.0';
const FO = 'urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0';
const XLINK = 'http://www.w3.org/1999/xlink';
const OFFICE = 'urn:oasis:names:tc:opendocument:xmlns:office:1.0';
const SVG = 'urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0';

const TAGS: Record<string, Record<string, string>> = {
  [TEXT]: {p: 'p', span: 'span', a: 'a', list: 'ul', 'list-item': 'li', 'list-header': 'li', 'line-break': 'br'},
  [TABLE]: {table: 'table', 'table-row': 'tr', 'table-cell': 'td'},
  [DRAW]: {image: 'img'},
};
const SKIPPED = new Set(['tracked-changes', 'sequence-decls', 'variable-decls', 'user-field-decls', 'note', 'annotation',
  'forms', 'table-columns', 'table-column', 'table-of-content-source', 'alphabetical-index-source',
  'illustration-index-source', 'bibliography-source', 'covered-table-cell']);

const imageType = (path: string) => {
  const ext = path.split('.').pop()!.toLowerCase();
  return `image/${ext === 'jpg' ? 'jpeg' : ext === 'svg' ? 'svg+xml' : ext}`;
};

/** Renders the text body of an OpenDocument text file (.odt) as HTML: headings, paragraphs, lists,
 * tables, links, images and the basic character and paragraph formatting of its styles. */
export async function renderOdt(bytes: Uint8Array, host: HTMLElement): Promise<void> {
  const zip = await JSZip.loadAsync(bytes);
  const parse = async (name: string) => {
    const xml = await zip.file(name)?.async('string');
    const doc = xml ? new DOMParser().parseFromString(xml, 'application/xml') : null;
    if (doc?.getElementsByTagName('parsererror').length)
      throw new Error(`${name} is not valid XML`);
    return doc;
  };
  const content = await parse('content.xml');
  if (!content)
    throw new Error('Not an OpenDocument file: content.xml is missing');
  const styles = new Map<string, Element>();
  for (const doc of [await parse('styles.xml'), content]) {
    for (const s of Array.from(doc?.getElementsByTagNameNS(STYLE, 'style') ?? []))
      styles.set(s.getAttributeNS(STYLE, 'name')!, s);
  }
  const images = new Map<string, string>();
  for (const frame of Array.from(content.getElementsByTagNameNS(DRAW, 'image'))) {
    const href = frame.getAttributeNS(XLINK, 'href');
    const file = href ? zip.file(href) : null;
    if (file && !images.has(href!))
      images.set(href!, `data:${imageType(href!)};base64,${await file.async('base64')}`);
  }

  const css = (styleName: string | null, seen = new Set<string>()): Record<string, string> => {
    const s = styleName ? styles.get(styleName) : undefined;
    if (!s || seen.has(styleName!))
      return {};
    seen.add(styleName!);
    const text = s.getElementsByTagNameNS(STYLE, 'text-properties')[0];
    const par = s.getElementsByTagNameNS(STYLE, 'paragraph-properties')[0];
    const cell = s.getElementsByTagNameNS(STYLE, 'table-cell-properties')[0];
    const fo = (e: Element | undefined, a: string) => e?.getAttributeNS(FO, a) || undefined;
    const underline = text?.getAttributeNS(STYLE, 'text-underline-style') || undefined;
    const align = fo(par, 'text-align');
    const own = {
      fontWeight: fo(text, 'font-weight'),
      fontStyle: fo(text, 'font-style'),
      fontSize: fo(text, 'font-size'),
      color: fo(text, 'color'),
      textDecoration: underline && (underline === 'none' ? 'none' : 'underline'),
      textAlign: align === 'start' ? 'left' : align === 'end' ? 'right' : align,
      marginLeft: fo(par, 'margin-left'),
      textIndent: fo(par, 'text-indent'),
      border: fo(cell, 'border'),
      borderTop: fo(cell, 'border-top'),
      borderBottom: fo(cell, 'border-bottom'),
      borderLeft: fo(cell, 'border-left'),
      borderRight: fo(cell, 'border-right'),
      padding: fo(cell, 'padding'),
    };
    return {...css(s.getAttributeNS(STYLE, 'parent-style-name'), seen),
      ...Object.fromEntries(Object.entries(own).filter(([_, v]) => v !== undefined))};
  };

  const convert = (node: Node, parent: HTMLElement): void => {
    if (node.nodeType === Node.TEXT_NODE) {
      parent.append(node.textContent!);
      return;
    }
    if (node.nodeType !== Node.ELEMENT_NODE)
      return;
    const e = node as Element;
    if (e.namespaceURI === SVG || SKIPPED.has(e.localName) ||
      (e.namespaceURI === DRAW && e.localName === 'image' && e.previousElementSibling?.localName === 'image'))
      return;
    if (e.namespaceURI === TEXT && (e.localName === 'tab' || e.localName === 's')) {
      parent.append('\u00a0'.repeat(e.localName === 'tab' ? 4 : Math.min(1000, +(e.getAttributeNS(TEXT, 'c') ?? 1) || 1)));
      return;
    }
    const tag = e.namespaceURI === TEXT && e.localName === 'h' ?
      `h${Math.min(6, Math.max(1, +(e.getAttributeNS(TEXT, 'outline-level') ?? 1) || 1))}` : TAGS[e.namespaceURI!]?.[e.localName];
    const el = tag ? document.createElement(tag) : parent;
    if (el !== parent) {
      Object.assign(el.style, css(e.getAttributeNS(TEXT, 'style-name') ?? e.getAttributeNS(TABLE, 'style-name')));
      const href = e.getAttributeNS(XLINK, 'href') ?? '';
      if (el instanceof HTMLAnchorElement && isSafeUrl(href))
        Object.assign(el, {href, target: '_blank', rel: 'noopener noreferrer'});
      if (el instanceof HTMLImageElement) {
        el.src = images.get(href) ?? '';
        el.style.width = e.parentElement?.getAttributeNS(SVG, 'width') ?? '';
        el.style.maxWidth = '100%';
      }
      if (el instanceof HTMLTableCellElement) {
        el.colSpan = +(e.getAttributeNS(TABLE, 'number-columns-spanned') ?? 1);
        el.rowSpan = +(e.getAttributeNS(TABLE, 'number-rows-spanned') ?? 1);
      }
      parent.append(el);
    }
    for (const child of Array.from(e.childNodes))
      convert(child, el);
  };

  const body = content.getElementsByTagNameNS(OFFICE, 'text')[0];
  for (const child of Array.from(body?.childNodes ?? []))
    convert(child, host);
}
