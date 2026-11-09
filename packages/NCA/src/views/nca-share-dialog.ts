/**
 * Share dialog — grant a group View access to a
 * study. The platform enforces it: the study row and, through it, its analyses and caches.
 *
 * View only: the platform does not yet honour a row Edit grant made from package
 * code. Snapshot files stay listable and readable by any user of the server until the platform
 * binds file access to row permissions (decision 9) — the dialog says so, so the scientist knows
 * what they rely on for IND-enabling data.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {getNcaRepository} from '../adapter';

/** Resolve a group typed as an id, a name or a friendly name. */
async function findGroup(text: string): Promise<{id: string; name: string} | null> {
  const groups = await grok.dapi.groups.list();
  const g = groups.find((x) => x.id === text || x.name === text || x.friendlyName === text);
  return g === undefined ? null : {id: g.id, name: g.friendlyName ?? g.name};
}

/** Open the share dialog for a study and grant the chosen group View access. */
export function openShareDialog(study: {id: string; label: string}, onDone: () => void): void {
  const groupInput = ui.input.string('Share with group', {
    value: '', nullable: true,
    tooltipText: 'The Datagrok group (name or id) that gets View access to this study.',
  });
  const note = ui.divText(
    'The group can view the study and its analyses; editing stays with you. Until a platform fix, ' +
    'snapshot files can be listed and read by any user of this server.',
    {classes: 'text-muted-foreground', style: {maxWidth: '440px', whiteSpace: 'normal'}});
  ui.dialog('Share study')
    .add(ui.divV([groupInput.root, note]))
    .onOK(() => {
      const text = (groupInput.value ?? '').trim();
      if (text.length === 0) return; // nothing to grant
      void (async () => {
        try {
          const group = await findGroup(text);
          if (group === null) {
            grok.shell.warning(`No group named '${text}'.`);
            return;
          }
          await getNcaRepository().share(study.id, group.id);
          grok.shell.info(`Shared '${study.label}' with '${group.name}' (view).`);
          onDone();
        } catch (e) {
          grok.shell.error(`Could not share '${study.label}': ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}
