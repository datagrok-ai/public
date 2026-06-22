import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface ExecError {
  blockIndex: number;
  error: string;
}

export async function executeSingleBlock(
  code: string, view: DG.ViewBase, blockIndex: number,
): Promise<{element: HTMLElement | null; error: ExecError | null}> {
  try {
    const t = view.type === DG.VIEW_TYPE.TABLE_VIEW ? (view as DG.TableView).dataFrame : undefined;
    const result = await new Function('grok', 'ui', 'DG', 'view', 't',
      'return (async () => {' + code + '})()',
    )(grok, ui, DG, view, t);
    return {element: result instanceof HTMLElement ? result : null, error: null};
  } catch (e: any) {
    grok.shell.error(`datagrok-exec error: ${e.message}`);
    return {element: null, error: {blockIndex, error: e?.message ?? String(e)}};
  }
}

export function buildViewContext(view: DG.ViewBase): string {
  if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
    const df = (view as DG.TableView).dataFrame;
    if (!df)
      return '';
    const cols = df.columns.toList().map((c) => `${c.name}(${c.type})`).join(', ');
    return `Table "${df.name}" (${df.rowCount} rows): ${cols}`;
  }
  if (view.type === 'ScriptView') {
    const scriptView = view as DG.ScriptView;
    const code = scriptView.code;
    if (code)
      return `ScriptView "${view.name}"\nCurrent script:\n\`\`\`\n${code}\n\`\`\``;
    return `ScriptView "${view.name}" (empty script)`;
  }
  return '';
}

interface DgEntityRef {
  type: 'file' | 'script' | 'query' | 'connection' | 'project' | 'space' | 'group' | 'user';
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
  case 'group':
    return ref.id ? grok.dapi.groups.find(ref.id) : null;
  case 'user':
    return ref.id ? grok.dapi.users.find(ref.id) : null;
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
