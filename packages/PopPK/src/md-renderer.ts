import * as ui from 'datagrok-api/ui';

type MdLinkHandler = (fileName: string) => void;
type DiagramRenderer = (src: string) => HTMLElement | null;

export interface RenderOptions {
  onMdLink?: MdLinkHandler;
  renderDiagram?: DiagramRenderer;
  glossary?: Map<string, string>;
}

export function renderMarkdown(text: string, options?: RenderOptions): HTMLElement {
  const opts = options ?? {};
  const lines = text.split('\n');
  const elements: HTMLElement[] = [];
  let i = 0;

  while (i < lines.length) {
    const line = lines[i];

    if (line.trim() === '') { i++; continue; }

    if (/^---+\s*$/.test(line.trim())) {
      elements.push(ui.element('hr'));
      i++;
      continue;
    }

    const headingMatch = line.match(/^(#{1,3})\s+(.*)/);
    if (headingMatch) {
      const level = headingMatch[1].length;
      const content = headingMatch[2];
      const heading = level === 1 ? ui.h1(content) : level === 2 ? ui.h2(content) : ui.h3(content);
      elements.push(heading);
      i++;
      continue;
    }

    if (line.includes('|') && i + 1 < lines.length && /^\s*\|?[\s:]*-+/.test(lines[i + 1])) {
      const tableEl = parseTable(lines, i, opts);
      elements.push(tableEl.element);
      i = tableEl.nextLine;
      continue;
    }

    if (line.trimStart().startsWith('>') || line.trimStart().startsWith('→')) {
      const quoteLines: string[] = [];
      while (i < lines.length && (lines[i].trimStart().startsWith('>') || lines[i].trimStart().startsWith('→'))) {
        quoteLines.push(lines[i].trimStart().replace(/^[>→]\s*/, ''));
        i++;
      }
      const quoteContent = renderInline(quoteLines.join(' '), opts);
      const blockquote = ui.div([quoteContent]);
      blockquote.style.borderLeft = '3px solid var(--grey-3, #ccc)';
      blockquote.style.paddingLeft = '12px';
      blockquote.style.margin = '8px 0';
      blockquote.style.color = 'var(--grey-5, #666)';
      blockquote.style.fontStyle = 'italic';
      elements.push(blockquote);
      continue;
    }

    if (/^\s*[-*]\s+/.test(line)) {
      const listItems: HTMLElement[] = [];
      while (i < lines.length && /^\s*[-*]\s+/.test(lines[i])) {
        const itemText = lines[i].replace(/^\s*[-*]\s+/, '');
        listItems.push(renderInline(itemText, opts));
        i++;
      }
      const ul = ui.element('ul');
      for (const item of listItems) {
        const li = ui.element('li');
        li.appendChild(item);
        ul.appendChild(li);
      }
      elements.push(ul);
      continue;
    }

    if (/^\s*\d+[.)]\s+/.test(line)) {
      const listItems: HTMLElement[] = [];
      while (i < lines.length && /^\s*\d+[.)]\s+/.test(lines[i])) {
        const itemText = lines[i].replace(/^\s*\d+[.)]\s+/, '');
        listItems.push(renderInline(itemText, opts));
        i++;
      }
      const ol = ui.element('ol');
      for (const item of listItems) {
        const li = ui.element('li');
        li.appendChild(item);
        ol.appendChild(li);
      }
      elements.push(ol);
      continue;
    }

    const imgMatch = line.match(/^!\[([^\]]*)\]\(([^)]+)\)/);
    if (imgMatch) {
      if (opts.renderDiagram) {
        const diagram = opts.renderDiagram(imgMatch[2]);
        if (diagram)
          elements.push(diagram);
      }
      i++;
      continue;
    }

    const paraLines: string[] = [];
    while (i < lines.length && lines[i].trim() !== '' &&
      !lines[i].match(/^#{1,3}\s/) && !lines[i].match(/^---+\s*$/) &&
      !(lines[i].includes('|') && i + 1 < lines.length && /^\s*\|?[\s:]*-+/.test(lines[i + 1])) &&
      !lines[i].trimStart().startsWith('>') && !lines[i].trimStart().startsWith('→') &&
      !/^\s*[-*]\s+/.test(lines[i]) && !/^\s*\d+[.)]\s+/.test(lines[i]) &&
      !lines[i].match(/^!\[/)) {
      paraLines.push(lines[i]);
      i++;
    }
    if (paraLines.length > 0) {
      const p = ui.div([renderInline(paraLines.join(' '), opts)]);
      p.style.margin = '6px 0';
      p.style.lineHeight = '1.5';
      elements.push(p);
    }
  }

  const container = ui.divV(elements);
  container.style.padding = '12px 20px';
  container.style.maxWidth = '960px';
  return container;
}

function parseTable(lines: string[], start: number, opts: RenderOptions):
  {element: HTMLElement; nextLine: number} {
  const parseCells = (line: string): string[] =>
    line.replace(/^\|/, '').replace(/\|$/, '').split('|').map((c) => c.trim());

  const headers = parseCells(lines[start]);
  let i = start + 2;
  const rows: string[][] = [];
  while (i < lines.length && lines[i].includes('|') && lines[i].trim() !== '') {
    rows.push(parseCells(lines[i]));
    i++;
  }

  const table = ui.element('table') as HTMLTableElement;
  table.style.borderCollapse = 'collapse';
  table.style.width = '100%';
  table.style.margin = '8px 0';
  table.style.fontSize = '13px';

  const cellStyle = (cell: HTMLElement, isHeader: boolean) => {
    cell.style.border = '1px solid var(--grey-2, #ddd)';
    cell.style.padding = '6px 10px';
    cell.style.textAlign = 'left';
    if (isHeader) {
      cell.style.fontWeight = 'bold';
      cell.style.backgroundColor = 'var(--grey-1, #f5f5f5)';
    }
  };

  const thead = ui.element('thead');
  const headerRow = ui.element('tr');
  for (const h of headers) {
    const th = ui.element('th');
    th.appendChild(renderInline(h, opts));
    cellStyle(th, true);
    headerRow.appendChild(th);
  }
  thead.appendChild(headerRow);
  table.appendChild(thead);

  const tbody = ui.element('tbody');
  for (const row of rows) {
    const tr = ui.element('tr');
    for (let c = 0; c < headers.length; c++) {
      const td = ui.element('td');
      td.appendChild(renderInline(row[c] ?? '', opts));
      cellStyle(td, false);
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);

  const wrapper = ui.div([table]);
  wrapper.style.overflowX = 'auto';
  return {element: wrapper, nextLine: i};
}

function renderTokens(tokens: Token[], opts: RenderOptions): HTMLElement {
  const span = ui.element('span');

  for (const part of tokens) {
    switch (part.type) {
    case 'bold-italic': {
      const bi = ui.element('strong');
      const biI = ui.element('em');
      if (part.children && hasNested(part.children))
        biI.appendChild(renderTokens(part.children, opts));
      else
        biI.textContent = part.text;
      bi.appendChild(biI);
      span.appendChild(bi);
      break;
    }
    case 'bold': {
      const b = ui.element('strong');
      if (part.children && hasNested(part.children))
        b.appendChild(renderTokens(part.children, opts));
      else
        b.textContent = part.text;
      span.appendChild(b);
      break;
    }
    case 'italic': {
      const em = ui.element('em');
      if (part.children && hasNested(part.children))
        em.appendChild(renderTokens(part.children, opts));
      else
        em.textContent = part.text;
      span.appendChild(em);
      break;
    }
    case 'code': {
      const code = ui.element('code');
      code.textContent = part.text;
      code.style.backgroundColor = 'var(--grey-1, #f0f0f0)';
      code.style.padding = '1px 4px';
      code.style.borderRadius = '3px';
      code.style.fontFamily = 'monospace';
      code.style.fontSize = '0.9em';
      span.appendChild(code);
      break;
    }
    case 'md-link': {
      const link = ui.link(part.text, () => opts.onMdLink?.(part.href!), part.href);
      span.appendChild(link);
      break;
    }
    case 'ext-link': {
      const a = ui.element('a') as HTMLAnchorElement;
      a.textContent = part.text;
      a.href = part.href!;
      a.target = '_blank';
      a.rel = 'noopener';
      a.classList.add('ui-link');
      span.appendChild(a);
      break;
    }
    default:
      renderTextWithGlossary(part.text, opts.glossary, span);
    }
  }

  return span;
}

function renderTextWithGlossary(text: string, glossaryMap: Map<string, string> | undefined, parent: HTMLElement): void {
  if (!glossaryMap || glossaryMap.size === 0) {
    parent.appendChild(document.createTextNode(text));
    return;
  }

  // Build regex from glossary keys, sorted longest-first to avoid partial matches
  const terms = Array.from(glossaryMap.keys()).sort((a, b) => b.length - a.length);
  const escaped = terms.map((t) => t.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'));
  const re = new RegExp(`\\b(${escaped.join('|')})\\b`, 'g');

  let last = 0;
  let m: RegExpExecArray | null;

  while ((m = re.exec(text)) !== null) {
    if (m.index > last)
      parent.appendChild(document.createTextNode(text.slice(last, m.index)));

    const term = m[1];
    const definition = glossaryMap.get(term)!;
    const termSpan = ui.element('span');
    termSpan.textContent = term;
    termSpan.style.borderBottom = '1px dotted var(--grey-4, #999)';
    termSpan.style.cursor = 'help';
    if (definition.includes('\n')) {
      const parts = definition.split('\n');
      const tooltipEl = ui.divV(parts.map((p, idx) => {
        const line = ui.element('div');
        line.textContent = p;
        if (idx === 0)
          line.style.cssText = 'font-weight:bold;margin-bottom:4px';
        else
          line.style.cssText = 'font-style:italic;opacity:0.9';
        return line;
      }));
      tooltipEl.style.maxWidth = '360px';
      ui.tooltip.bind(termSpan, () => tooltipEl);
    } else {
      ui.tooltip.bind(termSpan, definition);
    }
    parent.appendChild(termSpan);

    last = m.index + m[0].length;
  }

  if (last < text.length)
    parent.appendChild(document.createTextNode(text.slice(last)));
}

function hasNested(tokens: Token[]): boolean {
  return tokens.some((t) => t.type !== 'text');
}

function renderInline(text: string, opts: RenderOptions): HTMLElement {
  return renderTokens(tokenize(text), opts);
}

interface Token {
  type: 'text' | 'bold' | 'italic' | 'bold-italic' | 'code' | 'md-link' | 'ext-link';
  text: string;
  href?: string;
  children?: Token[];
}

function tokenize(text: string): Token[] {
  const tokens: Token[] = [];
  const re = /\[([^\]]+)\]\(([^)]+)\)|`([^`]+)`|\*{3}(.+?)\*{3}|\*{2}(.+?)\*{2}|(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)/g;
  let last = 0;
  let m: RegExpExecArray | null;

  while ((m = re.exec(text)) !== null) {
    if (m.index > last)
      tokens.push({type: 'text', text: text.slice(last, m.index)});

    if (m[1] != null) {
      const href = m[2];
      const linkType = href.endsWith('.md') ? 'md-link' : 'ext-link';
      tokens.push({type: linkType, text: m[1], href});
    } else if (m[3] != null) {
      tokens.push({type: 'code', text: m[3]});
    } else if (m[4] != null) {
      tokens.push({type: 'bold-italic', text: m[4], children: tokenize(m[4])});
    } else if (m[5] != null) {
      tokens.push({type: 'bold', text: m[5], children: tokenize(m[5])});
    } else if (m[6] != null) {
      tokens.push({type: 'italic', text: m[6], children: tokenize(m[6])});
    }
    last = m.index + m[0].length;
  }

  if (last < text.length)
    tokens.push({type: 'text', text: text.slice(last)});

  return tokens;
}
