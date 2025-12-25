/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {TodoView} from "./app/todo_app";
import '../css/styles.css';
export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


//name: todoApplication
//meta.role: app
//input: string path { optional: true; meta.url: true }
//output: view result
//meta.url: /todo
export function todoApplication(path?: string) {
  const parent = grok.functions.getCurrentCall();
  return new TodoView(parent);
}
