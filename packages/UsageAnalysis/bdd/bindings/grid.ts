/* The two elements the grid's context-menu feature reads that no platform kind names: the HELM
   web editor the Helm package opens from the menu's Current Value group. The editor is a headerless
   full-screen dialog, so it has no title to address it by; its tabs are plain buttons without the
   platform's tab-handle markup, and its HELM tab shows the notation in a contenteditable pane the
   editor marks with a test id. */
import {element} from '@datagrok-libraries/bdd';

element('HELM notation tab', {selector: '.hw-app__tab[data-pane="helm"]',
  description: 'the "HELM" tab of the open HELM web editor'});
element('HELM notation', {selector: '[data-testid="notation-pane-content"]',
  description: 'the HELM text of the open HELM web editor, shown by its HELM notation tab'});
