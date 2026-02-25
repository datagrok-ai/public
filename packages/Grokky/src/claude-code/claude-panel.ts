import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function executeDatagrokBlocks(content: string, view: DG.ViewBase): void {
  const re = /```datagrok-exec\n([\s\S]*?)```/g;
  let match: RegExpExecArray | null;
  while ((match = re.exec(content)) !== null) {
    const code = match[1];
    try {
      const t = view.type === DG.VIEW_TYPE.TABLE_VIEW ? (view as DG.TableView).dataFrame : undefined;
      new Function('grok', 'ui', 'DG', 'view', 't',
        'return (async () => {' + code + '})()',
      )(grok, ui, DG, view, t)
        .catch((e: any) => grok.shell.error(`datagrok-exec error: ${e.message}`));
    }
    catch (e: any) {
      grok.shell.error(`datagrok-exec error: ${e.message}`);
    }
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
  if (view.type === 'ScriptView')
    return `ScriptView "${view.name}"`;
  return '';
}
