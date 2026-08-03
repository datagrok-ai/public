/** The non-invasive entry point to the guide system: a small floating help
 *  button that opens a menu of tutorials and how-to answers. */

import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {GuideHost} from './guide-model';
import {GuideRunner} from './guide-runner';
import {TUTORIALS, QUESTIONS} from './guide-content';
import {setTid} from '../utils/test-ids';

/** A round, bottom-left floating "?" button. Visible but out of the way. */
export function createHelpButton(onClick: (ev: MouseEvent) => void): HTMLElement {
  const fab = ui.div([ui.iconFA('graduation-cap')], 'ff-help-fab');
  setTid(fab, 'help-button');
  ui.tooltip.bind(fab, 'Tutorials & help');
  fab.addEventListener('click', (ev) => onClick(ev));
  return fab;
}

/** Popup menu: Tutorials (guided, multi-step) and How do I…? (single answers). */
export function openGuideMenu(host: GuideHost, runner: GuideRunner, ev?: MouseEvent): void {
  const menu = DG.Menu.popup();
  const tutGroup = menu.group('Tutorials');
  for (const g of TUTORIALS)
    tutGroup.item(g.title, () => void runner.run(g, host));
  const howDoIGroup = menu.group('How do I…?');
  for (const g of QUESTIONS)
    howDoIGroup.item(g.title, () => void runner.run(g, host));
  menu.show(ev ? {causedBy: ev} : undefined);
}
