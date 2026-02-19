import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {MpoProfileInfo} from './utils';
import {MpoProfileCreateView} from './mpo-create-profile';
import {MpoProfileManager} from './mpo-profile-manager';

export class MpoProfileHandler extends DG.ObjectHandler<MpoProfileInfo> {
  get type(): string {
    return 'mpo-profile';
  }

  isApplicable(x: any): boolean {
    return x != null &&
      typeof x === 'object' &&
      typeof x.name === 'string' &&
      typeof x.properties === 'object' &&
      x.properties != null;
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

    const editBtn = ui.bigButton('Edit', () => MpoProfileHandler.edit(profile));
    const ribbon = ui.divH([editBtn], {style: {justifyContent: 'flex-end'}});
    ribbon.classList.add('d4-ribbon-item');

    const panel = ui.accordion();
    panel.addTitle(ui.label('MPO Profile'));
    panel.root.append(editor.root, ribbon);

    return panel.root;
  }

  constructor() {
    super();
    this.registerParamFunc('Edit', (p: MpoProfileInfo) => MpoProfileHandler.edit(p));
    this.registerParamFunc('Clone', (p: MpoProfileInfo) => MpoProfileHandler.clone(p));
    this.registerParamFunc('Delete', (p: MpoProfileInfo) => MpoProfileHandler.delete(p));
  }

  static edit(profile: MpoProfileInfo): void {
    const editable = structuredClone(profile);
    const view = new MpoProfileCreateView(editable, false, profile.fileName);
    grok.shell.v = grok.shell.addView(view.view);
  }

  static clone(profile: MpoProfileInfo): void {
    const {profile: clonedProfile, fileName} = MpoProfileManager.prepareClone(profile);
    const view = new MpoProfileCreateView(clonedProfile, false, fileName);
    grok.shell.v = grok.shell.addView(view.view);
  }

  static delete(profile: MpoProfileInfo): void {
    MpoProfileManager.confirmDelete(profile);
  }
}
