/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  InputForm, Viewer, DGBigButton, 
  DGButton, DGIconFA, 
  DGToggleInput, DGComboPopup, DGMarkdown,
  DGIcon
} from '@datagrok-libraries/webcomponents';
import {DockSpawnTsWebcomponent} from '@datagrok-libraries/webcomponents'

export const _package = new DG.Package();

let inited = false;

//tags: init
export function init() {
  if (inited) return;

  customElements.define('dg-viewer', Viewer);
  customElements.define('dg-input-form', InputForm);
  customElements.define('dg-button', DGButton, {extends: 'button'});
  customElements.define('dg-big-button', DGBigButton, {extends: 'button'});
  customElements.define('dg-icon-fa', DGIconFA);
  customElements.define('dg-icon', DGIcon);
  customElements.define('dg-toggle-input', DGToggleInput);
  customElements.define('dg-combo-popup', DGComboPopup);
  customElements.define('dock-spawn-ts', DockSpawnTsWebcomponent);
  customElements.define('dg-markdown', DGMarkdown);
  inited = true;
  console.log('webcomponents registered');
}
