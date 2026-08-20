import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {programGroupNames, T_ALIGNMENT, WORKSTREAMS} from '../domain/constants';
import {ensureProgram} from '../domain/security';
import {rejectPublication, updateCuration} from '../service/publication-service';
import {chips, compoundsChipInput, programCompoundCodes, tagsChipInput} from './registry-inputs';

/** Create/edit dialog for a program. Saves through {@link ensureProgram} — a plain
 * domain insert would skip the audience group provisioning, promote, and shares. */
export function showProgramDialog(existing?: {[key: string]: any}): DG.Dialog {
  const codeInput = ui.input.string('Code', {value: existing?.code ?? ''});
  codeInput.setTooltip('Unique program code; also seeds the audience group names');
  codeInput.enabled = existing == null;
  const nameInput = ui.input.string('Name', {value: existing?.name ?? ''});
  const descriptionInput = ui.input.textArea('Description', {value: existing?.description ?? ''});
  const externalRefInput = ui.input.string('External ref', {value: existing?.external_ref ?? '', nullable: true});
  externalRefInput.setTooltip('Link into the owning portfolio/DMS system');
  const stageInput = ui.input.string('Stage', {value: existing?.stage ?? '', nullable: true});
  const statusInput = ui.input.string('Status', {value: existing?.status ?? '', nullable: true});
  const indicationInput = ui.input.string('Indication', {value: existing?.indication ?? '', nullable: true});
  const groupsNote = ui.divText('', 'aa-program-groups-note');
  const groupLinks = ui.divV([], {style: {gap: '4px'}});
  if (existing != null) {
    groupLinks.appendChild(ui.divText('Audience groups:'));
    for (const id of [existing.viewers_group, existing.contributors_group, existing.approvers_group]) {
      if (id != null) {
        grok.dapi.groups.find(id)
          .then((group) => group != null && groupLinks.appendChild(ui.render(group)))
          .catch(() => {});
      }
    }
  }

  let dialog = ui.dialog(existing == null ? 'New program' : `Edit program ${existing.code}`);
  for (const input of [codeInput, nameInput, descriptionInput, externalRefInput,
    stageInput, statusInput, indicationInput])
    dialog = dialog.add(input);
  dialog = dialog.add(groupsNote).add(groupLinks)
    .onOK(async () => {
      const progress = DG.TaskBarProgressIndicator.create('Saving program…');
      try {
        const info = await ensureProgram({
          code: codeInput.value!.trim(),
          name: nameInput.value!.trim(),
          description: descriptionInput.value?.trim() || undefined,
          externalRef: externalRefInput.value?.trim() || undefined,
          stage: stageInput.value?.trim() || undefined,
          status: statusInput.value?.trim() || undefined,
          indication: indicationInput.value?.trim() || undefined,
        });
        grok.shell.info(`Program ${info.code} ${existing == null ? 'created' : 'updated'}`);
      } catch (e: any) {
        grok.shell.error(`Saving program failed: ${e?.message ?? e}`);
      } finally {
        progress.close();
      }
    })
    .show({center: true, width: 450});

  const okButton = dialog.getButton('OK');
  const refresh = () => {
    const code = codeInput.value?.trim();
    okButton.disabled = !code || !nameInput.value?.trim();
    const names = existing == null && code ? programGroupNames(code) : null;
    groupsNote.textContent = names == null ? '' :
      `Audience groups: ${names.viewers}, ${names.contributors}, ${names.approvers}`;
  };
  refresh();
  codeInput.onChanged.subscribe(refresh);
  nameInput.onChanged.subscribe(refresh);
  return dialog;
}

export function showRejectDialog(rowId: string): DG.Dialog {
  const reasonInput = ui.input.textArea('Reason', {value: ''});
  const dialog = ui.dialog('Reject publication')
    .add(reasonInput)
    .onOK(async () => {
      try {
        await rejectPublication(rowId, reasonInput.value!.trim());
        grok.shell.info('Publication rejected — reopen the view to refresh');
      } catch (e: any) {
        grok.shell.error(`Reject failed: ${e?.message ?? e}`);
      }
    })
    .show({center: true, width: 400});
  const okButton = dialog.getButton('OK');
  const refresh = () => okButton.disabled = !reasonInput.value?.trim();
  refresh();
  reasonInput.onChanged.subscribe(refresh);
  return dialog;
}

/** Edits the informational columns of a live version row; identity (name,
 * program, study) and the approval columns are deliberately absent. */
export async function showCurateDialog(rowId: string): Promise<DG.Dialog> {
  const [row] = await grok.dapi.domains.table(T_ALIGNMENT).query({
    filter: [{property: 'id', operator: '=', value: rowId}], expand: ['tags', 'compounds']});
  const linked = row?.program_id == null ? [] : await programCompoundCodes(row.program_id);
  const compoundsInput = compoundsChipInput(() => linked,
    (row?.compounds ?? []).map((c: any) => c.name));
  const tagsInput = tagsChipInput((row?.tags ?? []).map((t: any) => t.name));
  const workstreamInput = ui.input.choice<string | null>('Workstream',
    {items: [...WORKSTREAMS], value: row?.workstream ?? null, nullable: true});
  const pathInput = ui.input.string('Path', {value: row?.path ?? '', nullable: true});
  pathInput.setTooltip('Folder-like grouping, e.g. PK/exposure');
  const descriptionInput = ui.input.textArea('Description', {value: row?.description ?? ''});
  return ui.dialog('Curate publication')
    .add(ui.divV([compoundsInput, tagsInput, workstreamInput, pathInput, descriptionInput]))
    .onOK(async () => {
      try {
        await updateCuration(rowId, {
          path: pathInput.value?.trim() ?? '',
          tags: chips(tagsInput.value),
          compounds: chips(compoundsInput.value),
          description: descriptionInput.value?.trim() ?? '',
          workstream: workstreamInput.value ?? undefined,
        });
        grok.shell.info('Curation saved');
      } catch (e: any) {
        grok.shell.error(`Curation failed: ${e?.message ?? e}`);
      }
    })
    .show({center: true, width: 450});
}
