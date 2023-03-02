import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export function getMacromoleculeColumn(): DG.Column | any {
  const col = grok.shell.t.columns.bySemType(DG.SEMTYPE.MACROMOLECULE);
  if (col === null) {
    grok.shell.error('Current table does not contain macromolecules');
    return;
  }
  return col;
}

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}