import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {u2} from '@datagrok-libraries/utils/src/u2';
import {_package, PackageFunctions} from '../package';
import {MPO_TEMPLATE_PATH} from '../analysis/mpo';
import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';

export class InfoView {
  name: string;
  root: HTMLElement = ui.divV([]);

  constructor() {
    this.name = 'MPO Profiles';
    grok.shell.windows.showHelp = true;
  }

  onActivated(): void {
    grok.shell.windows.showHelp = true;
  }

  async init(): Promise<void> {
    ui.setUpdateIndicator(this.root, true);
    try {
      const appHeader = u2.appHeader({
        iconPath: _package.webRoot + '/images/icons/hit-triage-icon.png',
        learnMoreUrl: 'https://github.com/datagrok-ai/public/blob/master/packages/HitTriage/README_HT.md',
        description: '- Create and manage MPO profiles.\n'+
        '- Lorem ipsum.\n'+
        '- Lorem ipsum.\n'+
        '- Lorem ipsum.\n'+
        '- Lorem ipsum.\n ',
      });

      const manageProfilesHeader = ui.h1('Manage Profiles');
      const createNewProfileHeader = ui.h1('New Profile');
      const profilesTable = await this.getProfilesTable();

      const profilesDiv = ui.div([], {classes: 'ui-form'});
      await this.createProfileForm(profilesDiv);
      $(this.root).empty();
      this.root.appendChild(ui.divV([
        appHeader,
        manageProfilesHeader,
        profilesTable,
        createNewProfileHeader,
        profilesDiv,
      ]));
    } catch (e) {
      ui.setUpdateIndicator(this.root, false);
      throw e;
    } finally {
      ui.setUpdateIndicator(this.root, false);
    }
  }

  private async getProfilesTable() {
    const profileFiles = await grok.dapi.files.list(MPO_TEMPLATE_PATH);

    const profilesInfo = await Promise.all(profileFiles.map(async (file) => {
      const content =
        JSON.parse(await grok.dapi.files.readAsText(`${MPO_TEMPLATE_PATH}/${file.name}`)) as DesirabilityProfile;
      return {
        file: file,
        fileName: file.name,
        name: content.name || file.name.replace('.json', ''),
        description: content.description || '',
        properties: content.properties,
      };
    }));

    const deleteProfileIcon = (info: typeof profilesInfo[0]) => {
      const icon = ui.icons.delete(async () => {
        ui.dialog('Delete profile')
          .add(ui.divText(`Are you sure you want to delete profile "${info.name}"?`))
          .onOK(async () => {
            await grok.dapi.files.delete(`${MPO_TEMPLATE_PATH}/${info.fileName}`);
            await this.init();
          })
          .show();
      }, 'Delete profile');
      return icon;
    };

    const table = ui.table(profilesInfo, (info) => [
      ui.link(info.name, () => grok.shell.addPreview(PackageFunctions.mpoProfileEditor(info.file))),
      info.description,
      deleteProfileIcon(info),
    ], ['Name', 'Description', '']);

    table.style.color = 'var(--grey-5)';
    return table;
  }

  private async createProfileForm(containerDiv: HTMLElement) {
    $(containerDiv).empty();

    const methodChoice = ui.input.choice('Method', {
      items: ['Manual', 'Probabilistic'], value: 'Manual', nullable: false,
      tooltipText: 'To be added a nice tooltip text that explains both concepts',
    });

    const nameInput = ui.input.string('Profile Name', {tooltipText: 'to be added'});
    const descInput = ui.input.string('Description', {tooltipText: 'to be added'});

    // will need to add a validation because dataset is required for the pMPO approach (!important: labeled dataset)
    const datasetInput = ui.input.file('Dataset', {nullable: true});

    const createButton = ui.button('Create', async () => {
      const profileName = nameInput.value.trim();
      const description = descInput.value.trim();
      const method = methodChoice.value;
      const dataset = datasetInput.value;

      if (method === 'Manual')
        console.log(`Creating manual profile "${profileName}"`);
      else
        console.log(`Creating probabilistic profile "${profileName}" with dataset:`, dataset);
    });

    const form = ui.divV([
      methodChoice.root,
      nameInput.root,
      descInput.root,
      datasetInput.root,
      createButton,
    ], {style: {gap: '10px'}});

    containerDiv.appendChild(form);
  }
};
