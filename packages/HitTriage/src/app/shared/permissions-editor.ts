/* Edit Permissions: who may view a campaign and who may edit it, one row per user or group (u2 AccessGrid). The
   campaign keeps {view, edit} lists of group ids: edit implies view, and nobody listed means all users may view. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {AccessGrid, Dialog, computed, divV, signal, statusLine} from '@datagrok-libraries/u2';
import type {AccessRow} from '@datagrok-libraries/u2';
import {groupIcon, groupInput} from '@datagrok-libraries/u2/dg';
import {TriagePermissions} from '../types';

export const defaultPermissions: TriagePermissions = {
  view: [DG.Group.defaultGroupsIds['All users']],
  edit: [DG.Group.defaultGroupsIds['All users']],
};

export interface PermissionsEditorOptions {
  /** The campaign as the footer names it: "KHL-3 (BTK covalent follow-up)". */
  subject: string;
  /** Who edits it whatever the grid says: its author. */
  author?: string;
  /** Runs on OK, the dialog busy until it is done; a failure stays in the dialog. */
  save: (permissions: TriagePermissions) => unknown;
}

export class PermissionsEditor {
  static async show(permissions: TriagePermissions | undefined, options: PermissionsEditorOptions): Promise<void> {
    const current = permissions ?? defaultPermissions;
    const ids = [...new Set([...current.view, ...current.edit])];
    const groups = (await Promise.all(ids.map((id) => grok.dapi.groups.find(id))))
      .filter((g): g is DG.Group => g != null);
    let shown: AccessRow[] = groups.map((g) =>
      ({principal: g.id, can: {view: true, edit: current.edit.includes(g.id)}}));
    const known = new Map(groups.map((g) => [g.id, g]));
    const dialog = new Dialog('Edit Permissions', {name: 'edit-permissions'});
    dialog.build(() => {
      const grid: AccessGrid = new AccessGrid({
        name: 'access',
        fill: true,
        capabilities: [{name: 'view', label: 'Can View'}, {name: 'edit', label: 'Can Edit'}],
        principalLabel: 'User or group',
        principals: groups.map((g) => ({value: g.id, label: g.friendlyName})),
        icon: (id) => known.has(id) ? groupIcon(known.get(id)!) : null,
        value: shown,
        picker: (add) => groupInput({placeholder: '+ Add user or group…', onPick: (group, label) => {
          known.set(group.id, group);
          add({value: group.id, label});
        }}).root,
        onInput: (rows) => {
          const implied = PermissionsEditor.implied(shown, rows);
          shown = implied;
          if (implied !== rows)
            grid.value = implied;
        },
      });
      const rows = computed(() => grid.value ?? []);
      const editors = computed(() => rows.value.some((r) => r.can.edit));
      // what `stored` does with a list nobody can view: said only while the grid shows that
      const nobody = computed(() => !rows.value.some((r) => r.can.view || r.can.edit));
      dialog.add(divV([grid.root,
        statusLine(signal(null), computed(() => nobody.value ? 'With no one listed, all users can view.' :
          'Can Edit includes Can View.')),
        statusLine(signal(null), computed(() => editors.value ? `Saves the access of ${options.subject}` :
          `Only ${options.author ?? 'the author'} will be able to edit ${options.subject}`)),
      ]));
      dialog.onOK(() => options.save(PermissionsEditor.stored(grid.peek() ?? [])));
    });
    dialog.show({width: 440});
  }

  /** Checking Can Edit checks Can View; unchecking Can View unchecks Can Edit. Answers `next` when it holds. */
  static implied(previous: readonly AccessRow[], next: AccessRow[]): AccessRow[] {
    let changed = false;
    const rows = next.map((r) => {
      if (!r.can.edit || r.can.view)
        return r;
      changed = true;
      const viewWasOn = previous.find((p) => p.principal === r.principal)?.can.view === true;
      return {...r, can: viewWasOn ? {...r.can, edit: false} : {...r.can, view: true}};
    });
    return changed ? rows : next;
  }

  /** The grid as the campaign keeps it: editors can view too, and an empty view list is all users. */
  static stored(rows: readonly AccessRow[]): TriagePermissions {
    const edit = rows.filter((r) => r.can.edit).map((r) => r.principal);
    const view = rows.filter((r) => r.can.view || r.can.edit).map((r) => r.principal);
    return {edit, view: view.length > 0 ? view : [...defaultPermissions.view]};
  }
}
