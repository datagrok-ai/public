/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import {Subscription} from 'rxjs';
import {u2} from '@datagrok-libraries/utils/src/u2';

import {_package} from '../package';
import {MpoProfileInfo, updateMpoPath, MpoPathMode, MPO_PROFILES_NAME, MPO_PROFILE_CHANGED_EVENT, MPO_PROFILE_DELETED_EVENT} from './utils';
import {MpoProfileCreateView} from './mpo-create-profile';
import {mpoProfileStore} from './mpo-profile-store';
import {downloadProfile, uploadProfile} from './mpo-profile-actions';
import {MpoProfileHandler} from './mpo-profile-handler';
import {attachMpoProfilesAi} from '../ai-tools/mpo';

export class MpoProfilesView {
  name = MPO_PROFILES_NAME;
  root = ui.divV([]);
  view: DG.View;

  private tableContainer = ui.divV([]);
  private subs: Subscription[] = [];
  private previewedId: string | null = null;
  private descriptionObserver = new ResizeObserver((entries) => {
    for (const entry of entries) {
      const span = entry.target as HTMLElement;
      span.classList.toggle('chem-mpo-description-expandable',
        span.classList.contains('expanded') || span.scrollHeight > span.clientHeight);
    }
  });

  constructor() {
    this.view = DG.View.fromRoot(this.root);
    this.view.name = this.name;
    grok.shell.windows.showHelp = false;
    grok.shell.windows.showProperties = true;
    updateMpoPath(this.view, MpoPathMode.List);
    attachMpoProfilesAi(this);
  }

  async render(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      ui.empty(this.root);
      this.root.append(
        this.buildHeader(),
        ui.h1('Manage Profiles'),
        this.tableContainer,
        ui.divH(this.buildActionButtons()),
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
      await mpoProfileStore.ensureLoaded();
      this.rerenderTable();
    } finally {
      ui.setUpdateIndicator(this.tableContainer, false);
    }
  }

  private buildActionButtons(): HTMLElement[] {
    const createButton = ui.button(ui.h2('Create profile'), () => this.openCreateProfile());
    createButton.classList.add('chem-mpo-action-button');
    const importButton = ui.button(ui.h2('Upload'), () => uploadProfile());
    importButton.classList.add('chem-mpo-action-button');
    return [createButton, importButton];
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
    this.descriptionObserver.disconnect();
    ui.empty(this.tableContainer);

    if (mpoProfileStore.items.length === 0) {
      this.tableContainer.append(ui.h2('No MPO profiles yet'));
      this.tableContainer.append(ui.divText('An administrator can load the default profiles by running Chem:seedMpoProfiles.'));
      return;
    }

    this.tableContainer.append(this.buildProfilesTable());
  }

  private buildProfilesTable(): HTMLElement {
    const table = ui.table(
      mpoProfileStore.items,
      (profile) => [
        this.buildActionsButton(profile),
        this.buildProfileLink(profile),
        this.buildDescription(profile.description),
      ],
      ['', 'Name', 'Description'],
    );

    table.classList.add('chem-mpo-profiles-table');
    return table;
  }

  private buildProfileLink(profile: MpoProfileInfo): HTMLElement {
    const link = ui.link(profile.name, () => this.preview(profile));
    link.addEventListener('dblclick', () => MpoProfileHandler.edit(profile));
    return ui.bind(profile, link);
  }

  private buildActionsButton(profile: MpoProfileInfo): HTMLElement {
    const actionsButton = ui.button(
      '⋮',
      () => {
        const menu = ui.popupMenu();
        menu.item('Edit', () => MpoProfileHandler.edit(profile));
        menu.item('Clone', () => MpoProfileHandler.clone(profile));
        menu.item('Download', () => downloadProfile(profile));
        menu.separator().item('Delete', () => MpoProfileHandler.delete(profile));
        menu.show();
      },
      'Actions',
    );
    actionsButton.classList.add('chem-mpo-actions-button');
    return actionsButton;
  }

  private buildDescription(text: string): HTMLElement {
    const span = ui.divText(text, 'chem-mpo-description');
    span.onclick = () => {
      if (span.classList.contains('chem-mpo-description-expandable'))
        span.classList.toggle('expanded');
    };
    this.descriptionObserver.observe(span);
    return span;
  }

  preview(profile: MpoProfileInfo): void {
    this.previewedId = profile.id;
    grok.shell.windows.showContextPanel = true;
    grok.shell.o = profile;
  }

  openCreateProfile(name?: string): MpoProfileCreateView {
    // Clear stale context panel (e.g. a profile selected via ui.bind in the list)
    // until a dataset is loaded and MpoContextPanel takes over.
    grok.shell.o = null;
    const view = new MpoProfileCreateView();
    if (name)
      view.setProfileName(name);
    grok.shell.v = grok.shell.addPreview(view.tableView!);
    view.setupBreadcrumbs();
    return view;
  }

  private listenForChanges(): void {
    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_CHANGED_EVENT).subscribe(() => this.reloadProfiles()));
    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_DELETED_EVENT).subscribe((data) => {
      if (data?.id === this.previewedId) {
        grok.shell.o = null;
        this.previewedId = null;
      }
    }));
    this.subs.push(grok.events.onViewRemoving.subscribe((v) => {
      if (v.args.view.id === this.view.id)
        this.detach();
    }));
  }

  private detach(): void {
    this.descriptionObserver.disconnect();
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }

  async show() {
    await this.render();
    grok.shell.addPreview(this.view);
  }
}
