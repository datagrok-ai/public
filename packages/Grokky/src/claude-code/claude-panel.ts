import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export async function executeDatagrokBlocks(content: string, view: DG.ViewBase): Promise<HTMLElement[]> {
  const re = /```datagrok-exec\n([\s\S]*?)```/g;
  const results: HTMLElement[] = [];
  let match: RegExpExecArray | null;
  while ((match = re.exec(content)) !== null) {
    const code = match[1];
    try {
      const t = view.type === DG.VIEW_TYPE.TABLE_VIEW ? (view as DG.TableView).dataFrame : undefined;
      const result = await new Function('grok', 'ui', 'DG', 'view', 't',
        'return (async () => {' + code + '})()',
      )(grok, ui, DG, view, t);
      if (result instanceof HTMLElement)
        results.push(result);
    } catch (e: any) {
      grok.shell.error(`datagrok-exec error: ${e.message}`);
    }
  }
  return results;
}

export function buildViewContext(view: DG.ViewBase): string {
  if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
    const df = (view as DG.TableView).dataFrame;
    if (!df)
      return '';
    const cols = df.columns.toList().map((c) => `${c.name}(${c.type})`).join(', ');
    return `Table "${df.name}" (${df.rowCount} rows): ${cols}`;
  }
  if (view.type === 'ScriptView')
    return `ScriptView "${view.name}"`;
  return '';
}

interface DgEntityRef {
  type: 'file' | 'script' | 'query' | 'connection' | 'project' | 'space';
  name: string;
  id?: string;
  connector?: string;
  path?: string;
  isDirectory?: boolean;
}

async function fetchEntity(ref: DgEntityRef): Promise<any> {
  switch (ref.type) {
  case 'file':
    if (ref.isDirectory)
      return null;
    if (ref.connector && ref.name) {
      const files = await grok.dapi.files.list(`${ref.connector}/`, true, ref.name);
      return files.find((f) => f.fileName === ref.name) ?? null;
    }
    return null;
  case 'script':
    return ref.id ? grok.dapi.scripts.find(ref.id) : null;
  case 'query':
    return ref.id ? grok.dapi.queries.find(ref.id) : null;
  case 'connection':
    return ref.id ? grok.dapi.connections.find(ref.id) : null;
  case 'project':
    return ref.id ? grok.dapi.projects.find(ref.id) : null;
  case 'space':
    return ref.id ? grok.dapi.spaces.find(ref.id) : null;
  }
}

function renderEntityRef(ref: DgEntityRef): HTMLElement | null {
  const placeholder = ui.divText(ref.name);
  fetchEntity(ref).then((entity) => {
    if (entity && document.contains(placeholder))
      placeholder.replaceWith(ui.renderCard(entity));
  }).catch(() => {});
  return placeholder;
}

export function renderEntityBlocks(container: HTMLElement): void {
  for (const block of Array.from(container.querySelectorAll('code.language-datagrok-entities'))) {
    const pre = block.parentElement;
    if (!pre || pre.tagName !== 'PRE')
      continue;
    try {
      const refs: DgEntityRef[] = JSON.parse(block.textContent ?? '[]');
      if (!Array.isArray(refs) || refs.length === 0)
        continue;
      const cards = refs.map(renderEntityRef).filter((c): c is HTMLElement => c !== null);
      if (cards.length === 0)
        continue;
      const cardsContainer = ui.divV(cards, 'grokky-entity-cards');
      pre.replaceWith(cardsContainer);
    } catch (e) {
      console.warn('Failed to parse datagrok-entities block:', e);
    }
  }
}
