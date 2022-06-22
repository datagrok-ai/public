import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { c } from "./package";

export function sortObject(o) {
    return Object.keys(o).sort().reduce((r, k) => (r[k] = o[k], r), {});
}

export function updateDivInnerHTML(div: HTMLElement, content: any) {
    div.innerHTML = '';
    div.append(content);
  }

  export function addView(view: DG.ViewBase): DG.ViewBase {
    view.box = true;
    view.parentCall = c;
    view.path = '/' + view.name;
    grok.shell.addView(view);
    return view;
  }
