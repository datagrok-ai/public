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
import { MpoProfileCreateView } from './mpo-create-profile';

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
          ui.link('Create new profile', () => {
            const createdView = new MpoProfileCreateView();
            const v = grok.shell.addPreview(createdView.view);
            grok.shell.v = v;
            // grok.shell.v = view;
          }),
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
        ui.link(profile.name, () => {
          const editor = new MpoProfileEditor();
          editor.setProfile(profile);

          const panel = ui.accordion();
          const icon = ui.element('i');
          icon.className = 'grok-icon svg-icon svg-histogram';
          panel.addTitle(ui.span([icon, ui.label('MPO Profile')]));

          const editButton = ui.bigButton('Edit', () => {
            const editView = new MpoProfileCreateView(profile, false);
            const v = grok.shell.addView(editView.view);
            grok.shell.v = v;
            // const v = grok.shell.addView(PackageFunctions.mpoProfileEditor(profile.file));
            // grok.shell.v = v;
          });

          const buttonContainer = ui.divH([editButton], {style: {justifyContent: 'flex-end'}});
          buttonContainer.classList.add('d4-ribbon-item');

          panel.root.appendChild(editor.root);
          panel.root.appendChild(buttonContainer);

          grok.shell.o = panel.root;
        },
        ),
        (() => {
          const moreBtn = ui.button('...', () => handleBtnClick(), 'Click to ...');
          const handleBtnClick = () => {
            const menu = ui.popupMenu();
            menu.item('Delete', () => this.confirmAndDeleteProfile(profile));
            menu.show();
          };
          return moreBtn;
        })(),
        (() => {
          const descriptionDiv = ui.divText(profile.description, 'chem-mpo-description');
          descriptionDiv.onmouseenter = (e) => ui.tooltip.show(profile.description, e.x, e.y);
          descriptionDiv.onmouseleave = () => ui.tooltip.hide();
          return descriptionDiv;
        })(),
        // (() => {
        //   const icon = ui.icons.delete(
        //     () => this.confirmAndDeleteProfile(profile),
        //     'Delete profile',
        //   );
        //   icon.style.color = 'var(--blue-1)';
        //   return icon;
        // })(),
      ],
      ['Name', '', 'Description'/*, ''*/],
    );

    table.style.width = 'fit-content';
    table.style.color = 'var(--grey-5)';

    return table;
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
}
