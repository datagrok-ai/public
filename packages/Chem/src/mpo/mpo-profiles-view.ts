import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {u2} from '@datagrok-libraries/utils/src/u2';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {_package} from '../package';
import {MpoProfileInfo, loadMpoProfiles, deleteMpoProfile, createPathFromArr} from './utils';
import {MpoProfileCreateView} from './mpo-create-profile';

export class MpoProfilesView {
  name = 'MPO Profiles';
  root = ui.divV([]);
  view: DG.View;

  private tableContainer = ui.divV([]);
  private profiles: MpoProfileInfo[] = [];

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    // this.view.path = this.name;
    grok.shell.windows.showHelp = false;
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
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private async reloadProfiles(): Promise<void> {
    ui.setUpdateIndicator(this.tableContainer, true);
    try {
      this.profiles = await loadMpoProfiles();
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

    if (this.profiles.length === 0) {
      this.tableContainer.append(ui.h2('No MPO profiles yet'));
      return;
    }

    this.tableContainer.append(this.buildProfilesTable());
  }

  private buildProfilesTable(): HTMLElement {
    const table = ui.table(
      this.profiles,
      (profile) => [
        this.buildActionsButton(profile),
        ui.link(profile.name, () => this.openProfile(profile)),
        ui.divText(profile.description, 'chem-mpo-description'),
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

  private profileForEditing(profile: MpoProfileInfo): MpoProfileInfo {
    const copy = {...profile};
    delete (copy as any).file;
    return copy;
  }

  private openCloneProfile(profile: MpoProfileInfo): void {
    const clone = this.profileForEditing(profile);
    clone.name = `${profile.name} (Copy)`;

    const view = new MpoProfileCreateView(clone, false, profile.fileName);
    grok.shell.v = grok.shell.addView(view.view);
  }

  private openProfile(profile: MpoProfileInfo): void {
    const editor = new MpoProfileEditor();
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
    const editable = this.profileForEditing(profile);
    const view = new MpoProfileCreateView(editable, false, profile.fileName);
    grok.shell.v = grok.shell.addView(view.view);
  }

  private openCreateProfile(): void {
    const view = new MpoProfileCreateView();
    grok.shell.v = grok.shell.addPreview(view.tableView!);
  }

  private confirmDelete(profile: MpoProfileInfo): void {
    ui.dialog('Delete profile')
      .add(ui.divText(`Are you sure you want to delete profile "${profile.name}"?`))
      .onOK(async () => {
        try {
          await deleteMpoProfile(profile);
          this.profiles = this.profiles.filter((p) => p !== profile);
          this.rerenderTable();
        } catch (e) {
          grok.shell.error(`Failed to delete profile "${profile.name}": ${e instanceof Error ? e.message : e}`);
        }
      })
      .show();
  }

  async show() {
    await this.render();
    grok.shell.addPreview(this.view);
  }
}
