import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import {TriagePermissions} from '../types';
export const defaultPermissions: TriagePermissions = {
  view: [DG.Group.defaultGroupsIds['All users']],
  edit: [DG.Group.defaultGroupsIds['All users']],
};

export class PermissionsDialog {
  permissions: TriagePermissions;
  editGroupsInput: DG.InputBase<DG.Group[]>;
  viewGroupsInput: DG.InputBase<DG.Group[]>;

  constructor(permissions?: TriagePermissions) {
    this.permissions = permissions ?? defaultPermissions;
    // @ts-ignore TODO: fix typings
    this.editGroupsInput = ui.input.userGroups('Can Edit');
    this.editGroupsInput.root.style.marginBottom = '16px';
    // @ts-ignore TODO: fix typings
    this.viewGroupsInput = ui.input.userGroups('Can View');
  }

  private async getDialogContent(): Promise<HTMLElement> {
    const viewGroups: DG.Group[] = [];
    const editGroups: DG.Group[] = [];
    for (const groupId of this.permissions.view) {
      const group = await grok.dapi.groups.find(groupId);
      if (group)
        viewGroups.push(group);
    }
    for (const groupId of this.permissions.edit) {
      const group = await grok.dapi.groups.find(groupId);
      if (group)
        editGroups.push(group);
    }

    // values are set separately as there still are some bugs with input
    //@ts-ignore
    this.editGroupsInput.value = editGroups;
    //@ts-ignore
    this.viewGroupsInput.value = viewGroups;
    return ui.divV([
      this.viewGroupsInput.root,
      this.editGroupsInput.root,
    ]);
  }

  async show(onOk: (res: TriagePermissions) => void) {
    const content = await this.getDialogContent();
    const dialog = ui.dialog('Edit Permissions')
      .add(content)
      .onOK(() => {
        const editSet = new Set(this.editGroupsInput.value.map((group) => group.id));
        const viewSet = new Set(this.viewGroupsInput.value.map((group) => group.id));
        this.permissions.edit = Array.from(editSet);
        editSet.forEach((id) => viewSet.add(id)); // all who can edit can also view
        this.permissions.view = Array.from(viewSet);
        if (this.permissions.view.length === 0)
          this.permissions.view = defaultPermissions.view;
        onOk(this.permissions);
      });
    content.parentElement?.style && (content.parentElement.style.overflow = 'visible');
    dialog.show();
  }
}
