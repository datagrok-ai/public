import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package, PackageFunctions} from '../package';
import {MPO_TEMPLATE_PATH} from '../analysis/mpo';
import {DesirabilityProfile, PropertyDesirability} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

type MpoProfileInfo = {
  file: DG.FileInfo;
  fileName: string;
  name: string;
  description: string;
  properties: {[key: string]: PropertyDesirability};
};

export class MpoProfilesView {
  readonly name = 'MPO Profiles';
  readonly root: HTMLElement = ui.divV([]);

  constructor() {
    grok.shell.windows.showHelp = true;
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
  }

  async render(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      $(this.root).empty();
      this.root.appendChild(
        ui.divV([
          this.buildAppHeader(),
          ui.h1('Manage Profiles'),
          await this.buildProfilesTable(),
          ui.h1('New Profile'),
          this.buildCreateProfileFormContainer(),
        ]),
      );
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private buildAppHeader(): HTMLElement {
    return u2.appHeader({
      iconPath: `${_package.webRoot}/images/icons/hit-triage-icon.png`,
      learnMoreUrl:
        'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
      description:
        '- Create and manage MPO profiles.\n' +
        '- Lorem ipsum.\n' +
        '- Lorem ipsum.\n' +
        '- Lorem ipsum.\n' +
        '- Lorem ipsum.\n',
    });
  }

  private async buildProfilesTable(): Promise<HTMLElement> {
    const profiles = await this.loadProfiles();

    const table = ui.table(
      profiles,
      (profile) => [
        ui.link(profile.name, () =>
          grok.shell.addPreview(
            PackageFunctions.mpoProfileEditor(profile.file),
          ),
        ),
        profile.description,
        ui.icons.delete(
          () => this.confirmAndDeleteProfile(profile),
          'Delete profile',
        ),
      ],
      ['Name', 'Description', ''],
    );

    table.style.color = 'var(--grey-5)';
    return table;
  }

  private buildCreateProfileFormContainer(): HTMLElement {
    const container = ui.div([], {classes: 'ui-form'});
    container.appendChild(this.buildCreateProfileForm());
    return container;
  }

  private buildCreateProfileForm(): HTMLElement {
    const methodInput = ui.input.choice('Method', {
      items: ['Manual', 'Probabilistic'],
      value: 'Manual',
      nullable: false,
      tooltipText: 'Explanation of Manual vs Probabilistic MPO',
    });

    const nameInput = ui.input.string('Profile Name', {tooltipText: 'Name of the MPO profile'});

    const descriptionInput = ui.input.string('Description', {tooltipText: 'Optional description'});

    const datasetInput = ui.input.table('Dataset', {nullable: true});

    const createButton = ui.button('Create', async () => {
      if (!datasetInput.value) {
        ui.dialog('Error').add(ui.divText('Dataset is required')).show();
        return;
      }

      const df = datasetInput.value;
      const newProfile: DesirabilityProfile = {
        name: nameInput.value.trim() || 'New Profile',
        description: descriptionInput.value.trim(),
        properties: {},
      };

      Array.from(df.columns.numerical).forEach((col: DG.Column) => {
        newProfile.properties[col.name] = {
          weight: 0.5,
          min: col.min,
          max: col.max,
          line: [],
        };
      });

      this.openNewProfileView(newProfile, df);
    });

    return ui.divV(
      [
        methodInput.root,
        nameInput.root,
        descriptionInput.root,
        datasetInput.root,
        createButton,
      ],
      {style: {gap: '10px'}},
    );
  }

  private async loadProfiles(): Promise<MpoProfileInfo[]> {
    const files = await grok.dapi.files.list(MPO_TEMPLATE_PATH);

    return Promise.all(
      files.map(async (file) => {
        const content = JSON.parse(
          await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`),
        ) as DesirabilityProfile;

        return {
          file,
          fileName: file.name,
          name: content.name ?? file.name.replace('.json', ''),
          description: content.description ?? '',
          properties: content.properties,
        };
      }),
    );
  }

  private async confirmAndDeleteProfile(profile: MpoProfileInfo): Promise<void> {
    ui.dialog('Delete profile')
      .add(ui.divText(`Are you sure you want to delete profile "${profile.name}"?`))
      .onOK(async () => {
        await grok.dapi.files.delete(`${MPO_TEMPLATE_PATH}/${profile.fileName}`);
        await this.render();
      })
      .show();
  }

  private openNewProfileView(profile: DesirabilityProfile, df: DG.DataFrame) {
    const view = DG.View.create();
    view.name = profile.name;

    const saveButton = ui.bigButton('SAVE', () => saveProfile());
    saveButton.style.display = 'none';
    view.setRibbonPanels([[saveButton]]);

    const editor = new MpoProfileEditor(df);
    editor.setProfile(profile);

    editor.onChanged.subscribe(() => saveButton.style.display = 'initial');

    const saveProfile = () => {
      grok.dapi.files.writeAsText(`${profile.name}.json`, JSON.stringify(profile))
        .then(() => {
          ui.dialog('Saved').add(ui.divText(`Profile saved as ${profile.name}.json`)).show();
          saveButton.style.display = 'none';
          this.render();
        })
        .catch((err) => ui.dialog('Error').add(ui.divText(err.message)).show());
    };

    view.append(editor.root);
    grok.shell.addView(view);
  }
}
