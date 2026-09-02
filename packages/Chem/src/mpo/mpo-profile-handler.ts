import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {isDesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileInfo} from './utils';
import {MpoProfileCreateView} from './mpo-create-profile';
import {confirmDeleteProfile, downloadProfile, prepareProfileClone} from './mpo-profile-actions';

export class MpoProfileHandler extends DG.ObjectHandler<MpoProfileInfo> {
  get type(): string {
    return 'mpo-profile';
  }

  isApplicable(x: any): boolean {
    return isDesirabilityProfile(x);
  }

  renderTooltip(profile: MpoProfileInfo): HTMLElement {
    return ui.divV([
      ui.h3(profile.name),
      ui.divText(profile.description || 'No description'),
    ]);
  }

  renderProperties(profile: MpoProfileInfo): HTMLElement {
    const editor = new MpoProfileEditor(undefined, false, true);
    editor.setProfile(profile);
    editor.root.style.pointerEvents = 'none';

    const ribbon = ui.divH([ui.bigButton('Edit', () => MpoProfileHandler.edit(profile))],
      {style: {justifyContent: 'flex-end'}});
    ribbon.classList.add('d4-ribbon-item');

    const panel = ui.accordion();
    panel.addTitle(ui.label('MPO Profile'));
    panel.root.append(editor.root, ribbon);

    return panel.root;
  }

  constructor() {
    super();
    this.registerParamFunc('Edit', MpoProfileHandler.edit);
    this.registerParamFunc('Clone', MpoProfileHandler.clone);
    this.registerParamFunc('Download', downloadProfile);
    this.registerParamFunc('Delete', MpoProfileHandler.delete);
  }

  static edit(profile: MpoProfileInfo): void {
    const view = new MpoProfileCreateView(profile, false);
    grok.shell.v = grok.shell.addView(view.view);
    view.setupBreadcrumbs();
  }

  static clone(profile: MpoProfileInfo): void {
    const view = new MpoProfileCreateView(prepareProfileClone(profile), false);
    grok.shell.v = grok.shell.addView(view.view);
    view.setupBreadcrumbs();
  }

  static delete(profile: MpoProfileInfo): void {
    confirmDeleteProfile(profile);
  }
}
