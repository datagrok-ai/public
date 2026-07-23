import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export interface ExecError {
  blockIndex: number;
  error: string;
}

export async function executeSingleBlock(
  code: string, view: DG.ViewBase, blockIndex: number,
): Promise<{element: HTMLElement | null; value: any; error: ExecError | null}> {
  try {
    // Prefer the live current view: Claude may have opened a new view earlier in the same turn.
    const liveView = grok.shell.v ?? view;
    const t = liveView?.type === DG.VIEW_TYPE.TABLE_VIEW ? (liveView as DG.TableView).dataFrame : undefined;
    const result = await new Function('grok', 'ui', 'DG', 'view', 't',
      'return (async () => {' + code + '})()',
    )(grok, ui, DG, liveView, t);
    const element = result instanceof HTMLElement ? result : null;
    return {element, value: element ? undefined : result, error: null};
  } catch (e: any) {
    return {element: null, value: undefined, error: {blockIndex, error: e?.message ?? String(e)}};
  }
}

const VERIFY_TIMEOUT_MS = 30000;

export interface VerificationResult {
  passed: boolean;
  observed: any;
  error: string | null;
}

export async function runVerification(
  assertion: string, view: DG.ViewBase,
): Promise<VerificationResult> {
  const liveView = grok.shell.v ?? view;
  const timeout = new Promise<null>((resolve) => setTimeout(() => resolve(null), VERIFY_TIMEOUT_MS));
  const res = await Promise.race([executeSingleBlock(assertion, liveView, 0), timeout]);
  if (res === null)
    return {passed: false, observed: undefined, error: `Verification timed out after ${VERIFY_TIMEOUT_MS / 1000}s`};
  if (res.element)
    return {passed: false, observed: undefined, error: 'Assertion must return the observed value, not a DOM element'};
  return {passed: res.error == null && !!res.value, observed: res.value, error: res.error?.error ?? null};
}

function describeTable(df: DG.DataFrame): string {
  const cols = df.columns.toList().map((c) => `${c.name}(${c.type})`).join(', ');
  return `Table "${df.name}" (${df.rowCount} rows): ${cols}`;
}

export function buildViewContext(view: DG.ViewBase): string {
  if (view.type === DG.VIEW_TYPE.TABLE_VIEW) {
    const df = (view as DG.TableView).dataFrame;
    return df ? describeTable(df) : '';
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

/** Snapshot of everything the user has open: current view (detailed), other views, and workspace tables. */
export function buildWorkspaceContext(): string {
  const lines: string[] = [];
  const current = grok.shell.v;
  if (current) {
    lines.push(`Current view: "${current.name}" (${current.type})`);
    const ai = (current as any).aiDescription;
    if (ai)
      lines.push(`About this view: ${ai}`);
    const details = buildViewContext(current);
    if (details)
      lines.push(details);
  }
  const others = Array.from(grok.shell.views).filter((v) => v !== current);
  if (others.length > 0) {
    lines.push('Other open views: ' + others.map((v) => {
      const df = v.type === DG.VIEW_TYPE.TABLE_VIEW ? (v as DG.TableView).dataFrame : null;
      return `"${v.name}" (${v.type}${df ? `, table "${df.name}"` : ''})`;
    }).join(', '));
  }
  const currentTable = current?.type === DG.VIEW_TYPE.TABLE_VIEW ? (current as DG.TableView).dataFrame : null;
  for (const t of grok.shell.tables.filter((t) => t !== currentTable))
    lines.push(describeTable(t));
  return lines.length > 0 ? 'Workspace state (live, changes as the user navigates):\n' + lines.join('\n') : '';
}

export interface DgEntityRef {
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

export function renderEntityRefList(refs: DgEntityRef[]): HTMLElement {
  const cards = refs.map(renderEntityRef).filter((c): c is HTMLElement => c !== null);
  return ui.divV(cards, 'grokky-entity-cards');
}

