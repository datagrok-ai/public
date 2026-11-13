/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  InputForm, Viewer, DGBigButton,
  DGButton, DGIconFA,
  DGToggleInput, DGComboPopup, DGMarkdown,
  ValidationIcon,
  DGIconImage,
} from '@datagrok-libraries/webcomponents';
import {DockSpawnTsWebcomponent} from '@datagrok-libraries/dock-spawn-dg';
import './styles.css';
export * from './package.g';
export const _package = new DG.Package();

const ELEMENTS: [string, CustomElementConstructor, ElementDefinitionOptions?][] = [
  ['dg-viewer', Viewer],
  ['dg-input-form', InputForm],
  ['dg-button', DGButton, {extends: 'button'}],
  ['dg-big-button', DGBigButton, {extends: 'button'}],
  ['dg-icon-fa', DGIconFA],
  ['dg-icon-image', DGIconImage],
  ['dg-toggle-input', DGToggleInput],
  ['dg-combo-popup', DGComboPopup],
  ['dock-spawn-ts', DockSpawnTsWebcomponent],
  ['dg-markdown', DGMarkdown],
  ['dg-validation-icon', ValidationIcon],
];

export class PackageFunctions {
  @grok.decorators.init()
  static async init() {
    for (const [name, element, options] of ELEMENTS) {
      if (!customElements.get(name))
        customElements.define(name, element, options);
    }
  }
}
