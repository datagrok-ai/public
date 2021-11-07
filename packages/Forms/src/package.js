/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok'
import * as ui from 'datagrok-api/ui'
import * as DG from 'datagrok-api/dg'

/* 3rd party imports */

/* Project imports */
import Editor from './Editor'

/* constants */
export const _package = new DG.Package()

//name: FormsApp
//tags: app
export function FormsApp() {
  const view = grok.shell.newView('Forms')
  const rootEl = view.root
  rootEl.setAttribute('id', 'gjs-root')
  rootEl.style.backgroundColor = '#2a2a2a'

  const windows = grok.shell.windows;
  windows.showToolbox = false;
  windows.showProperties = true;
  windows.showHelp = false;

  const editor = new Editor({
    view
  })
}
