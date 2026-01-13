import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package, PackageFunctions} from '../package';
import {MPO_TEMPLATE_PATH} from '../analysis/mpo';
import {DesirabilityProfile, PropertyDesirability} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MpoContextPanel} from './mpo-context-panel';

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

  private profilesTableContainer: HTMLElement = ui.divV([]);

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
          this.profilesTableContainer,
          ui.h1('New Profile'),
          this.createProfileForm(),
        ]),
      );

      ui.setUpdateIndicator(this.profilesTableContainer, true);
      const table = await this.buildProfilesTable();
      $(this.profilesTableContainer).empty().append(table);
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private buildAppHeader(): HTMLElement {
    return u2.appHeader({
      iconPath: `${_package.webRoot}/images/mpo.png`,
      learnMoreUrl:
        '',
      description:
        '- Create and manage MPO profiles.\n' +
        '- Build profiles manually or train from labeled data.\n' +
        '- Initialize profiles from scratch or from an existing dataset.\n' +
        '- Full lifecycle support: create, edit, clone, and delete profiles..\n',
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
        ui.divText(profile.description, 'chem-mpo-description'),
        (() => {
          const icon = ui.icons.delete(
            () => this.confirmAndDeleteProfile(profile),
            'Delete profile',
          );
          icon.style.color = 'var(--blue-1)';
          return icon;
        })(),
      ],
      ['Name', 'Description', ''],
    );
    table.style.width = 'fit-content';
    table.style.color = 'var(--grey-5)';

    return table;
  }

  private createProfileForm(): HTMLElement {
    const methodInput = ui.input.choice('Method', {
      items: ['Manual', 'Probabilistic'],
      value: 'Manual',
      nullable: false,
      tooltipText: 'Explanation of Manual vs Probabilistic MPO',
      onValueChanged: (value) => {
        const isProbabilistic = value === 'Probabilistic';
        nameInput.classList.toggle('chem-mpo-d-none', isProbabilistic);
        descriptionInput.classList.toggle('chem-mpo-d-none', isProbabilistic);
        datasetInput.nullable = !isProbabilistic;
      },
    });

    const nameInput = ui.input.string('Profile Name', {tooltipText: 'Name of the MPO profile', nullable: true});
    const descriptionInput = ui.input.string('Description', {tooltipText: 'Optional description', nullable: true});
    const datasetInput = ui.input.table('Dataset', {nullable: true});

    const createButton = ui.bigButton('Create', async () => {
      const df = datasetInput.value!;
      if (methodInput.value === 'Probabilistic') {
        this.trainPMPOProfile(df);
        return;
      }
      this.manualMPOProfile(df, nameInput.value, descriptionInput.value);
    });

    const buttonContainer = ui.divH([createButton], {style: {justifyContent: 'flex-end'}});
    buttonContainer.classList.add('d4-ribbon-item');

    const form = ui.form([
      methodInput,
      nameInput,
      descriptionInput,
      datasetInput,
    ]);

    return ui.divV([form, buttonContainer], {style: {gap: '10px'}});
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
        try {
          await grok.dapi.files.delete(`${MPO_TEMPLATE_PATH}/${profile.fileName}`);

          // Yeah, maybe it is better to remove from the table directly, than re-rendering all
          const rows = Array.from(this.profilesTableContainer.querySelectorAll('tr'));
          for (const row of rows) {
            if (row.textContent?.includes(profile.name)) {
              row.remove();
              break;
            }
          }
        } catch (e) {
          grok.shell.error(`Failed to delete profile "${profile.name}": ${e instanceof Error ? e.message : e}`);
        } finally {
        }
      })
      .show();
  }


  private async trainPMPOProfile(df: DG.DataFrame): Promise<void> {
    grok.shell.addTableView(df);
    try {
      await grok.functions.call('EDA:trainPmpo', {});
    } catch (err: any) {
      console.log(err);
    }
  }

  private async manualMPOProfile(
    df: DG.DataFrame | null, profileName: string | null, profileDescription: string | null,
  ): Promise<void> {
    const newProfile: DesirabilityProfile = {
      name: profileName ?? '',
      description: profileDescription ?? '',
      properties: {},
    };

    if (df) {
      const numericalColumns = Array.from(df.columns.numerical).slice(0, 3);
      numericalColumns.forEach((col: DG.Column) => {
        newProfile.properties[col.name] = {
          weight: 0.5,
          min: col.min,
          max: col.max,
          line: [],
        };
      });
      this.openNewProfileView(newProfile, df);
    } else {
      const view = DG.View.create();
      view.name = profileName ?? 'New MPO Profile';
      for (let i = 1; i <= 3; i++)
        newProfile.properties[`Property ${i}`] = {weight: 0.5, line: []};
      const editor = new MpoProfileEditor(undefined, true);
      editor.setProfile(newProfile);
      const scrollableEditor = ui.divV([editor.root], {style: {overflow: 'auto', height: '100%'}});
      view.append(scrollableEditor);
      grok.shell.addView(view);
    }
  }

  private openNewProfileView(profile: DesirabilityProfile, df: DG.DataFrame) {
    const view = grok.shell.addTableView(df);

    const saveButton = ui.bigButton('SAVE', () => saveProfile());
    saveButton.style.display = 'none';
    view.setRibbonPanels([[saveButton]]);

    const editor = new MpoProfileEditor(df, true);
    editor.setProfile(profile);

    editor.onChanged.subscribe(() => saveButton.style.display = 'initial');

    const saveProfile = () => {
      grok.dapi.files.writeAsText(`${profile.name}.json`, JSON.stringify(profile))
        .then(() => {
          grok.shell.info(`Profile "${profile.name}" saved successfully.`);
          saveButton.style.display = 'none';
          this.render();
        })
        .catch((err) => grok.shell.error(`Failed to save profile "${profile.name}": ${err}`));
    };
    const scrollableEditor = ui.divV([editor.root], {style: {overflow: 'auto', height: '100%'}});
    view.dockManager.dock(scrollableEditor, DG.DOCK_TYPE.DOWN, null, 'MPO Profile Editor', 0.6);

    const mpoContextPanel = new MpoContextPanel(df);
    grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      if (profile && mpoContextPanel) {
        await mpoContextPanel.render(
          profile,
          editor.columnMapping,
          'Average',
        );
      }
    });
    grok.shell.o = mpoContextPanel.root;
  }
}
