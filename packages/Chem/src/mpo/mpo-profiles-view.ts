/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {_package} from '../package';
import {MpoProfileInfo, updateMpoPath, MpoPathMode, MPO_PROFILE_CHANGED_EVENT} from './utils';
import {MpoProfileCreateView} from './mpo-create-profile';
import {MpoProfileManager} from './mpo-profile-manager';

export class MpoProfilesView {
  name = 'MPO Profiles';
  root = ui.divV([]);
  view: DG.View;

  private tableContainer = ui.divV([]);
  private subs: Subscription[] = [];

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    grok.shell.windows.showHelp = false;
    updateMpoPath(this.view, MpoPathMode.List);
  }

  async render(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    const createButton = ui.button(ui.h2('Create profile'), () => this.openCreateProfile());
    createButton.classList.add('chem-mpo-create-profile-button');
    try {
      ui.empty(this.root);
      this.root.append(
        this.buildHeader(),
        ui.h1('Manage Profiles'),
        this.tableContainer,
        createButton,
      );

      await this.reloadProfiles();
      this.listenForChanges();
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private async reloadProfiles(): Promise<void> {
    ui.setUpdateIndicator(this.tableContainer, true);
    try {
      await MpoProfileManager.ensureLoaded();
      this.rerenderTable();
    } finally {
      ui.setUpdateIndicator(this.tableContainer, false);
    }
  }

  private buildHeader(): HTMLElement {
    return u2.appHeader({
      iconPath: `${_package.webRoot}/images/mpo.png`,
      description:
        '- Create and manage MPO profiles.\n' +
        '- Build profiles manually or train from labeled data.\n' +
        '- Initialize profiles from scratch or from an existing dataset.\n' +
        '- Full lifecycle support: create, edit, clone, and delete profiles.\n',
    });
  }

  private rerenderTable(): void {
    ui.empty(this.tableContainer);

    if (MpoProfileManager.items.length === 0) {
      this.tableContainer.append(ui.h2('No MPO profiles yet'));
      return;
    }

    this.tableContainer.append(this.buildProfilesTable());
  }

  private buildProfilesTable(): HTMLElement {
    const table = ui.table(
      MpoProfileManager.items,
      (profile) => [
        this.buildActionsButton(profile),
        ui.link(profile.name, () => this.openProfile(profile)),
        this.buildDescription(profile.description),
      ],
      ['', 'Name', 'Description'],
    );

    table.classList.add('chem-mpo-profiles-table');
    return table;
  }

  private buildActionsButton(profile: MpoProfileInfo): HTMLElement {
    const actionsButton = ui.button(
      '⋮',
      () => {
        ui.popupMenu()
          .item('Edit', () => this.openEditProfile(profile))
          .item('Clone', () => this.openCloneProfile(profile))
          .item('Delete', () => this.confirmDelete(profile))
          .show();
      },
      'Actions',
    );
    actionsButton.classList.add('chem-mpo-actions-button');
    return actionsButton;
  }

  private buildDescription(text: string): HTMLElement {
    const span = ui.divText(text, 'chem-mpo-description');
    span.addEventListener('click', () => span.classList.toggle('expanded'));
    return span;
  }

  private openCloneProfile(profile: MpoProfileInfo): void {
    MpoProfileManager.clone(profile);
  }

  private openProfile(profile: MpoProfileInfo): void {
    const editor = new MpoProfileEditor(undefined, false, true);
    editor.setProfile(profile);
    editor.root.style.pointerEvents = 'none';

    const panel = ui.accordion();
    panel.addTitle(ui.label('MPO Profile'));
    panel.root.append(editor.root, this.buildEditRibbon(profile));

    grok.shell.o = panel.root;
  }

  private buildEditRibbon(profile: MpoProfileInfo): HTMLElement {
    const editBtn = ui.bigButton('Edit', () => this.openEditProfile(profile));

    const container = ui.divH([editBtn], {style: {justifyContent: 'flex-end'}});
    container.classList.add('d4-ribbon-item');
    return container;
  }

  private openEditProfile(profile: MpoProfileInfo): void {
    const editable = structuredClone(profile);
    const view = new MpoProfileCreateView(editable, false, profile.fileName);
    grok.shell.v = grok.shell.addView(view.view);
  }

  private openCreateProfile(): void {
    const view = new MpoProfileCreateView();
    grok.shell.v = grok.shell.addPreview(view.tableView!);
  }

  private confirmDelete(profile: MpoProfileInfo): void {
    MpoProfileManager.confirmDelete(profile);
  }

  private listenForChanges(): void {
    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_CHANGED_EVENT).subscribe(() => this.reloadProfiles()));
    this.subs.push(grok.events.onViewRemoving.subscribe((v) => {
      if (v.args.view.id === this.view.id)
        this.detach();
    }));
  }

  private detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }

  async show() {
    await this.render();
    grok.shell.addPreview(this.view);
  }
}
