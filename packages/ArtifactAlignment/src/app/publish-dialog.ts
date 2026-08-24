import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {T_ALIGNMENT, T_PROGRAM, T_STUDY, WORKSTREAMS} from '../domain/constants';
import {publishWorkflowRun} from '../service/publication-service';
import {chips, compoundsChipInput, programCompoundCodes, tagsChipInput} from './registry-inputs';

interface ProgramChoice {
  id: string;
  code: string;
  name: string;
}

/** The "Publish to program" dialog over a saved or live Compute2 run — doc § Publish
 * & republish flow. */
export async function showPublishDialog(source: string | DG.FuncCall, defaultName?: string):
  Promise<DG.Dialog | null> {
  const programs: ProgramChoice[] =
    await grok.dapi.domains.table(T_PROGRAM).query({sort: 'code', limit: 1000});
  if (programs.length === 0) {
    grok.shell.warning('No programs are visible to you — ask an administrator to provision one');
    return null;
  }

  const nameInput = ui.input.string('Name', {value: defaultName ?? '', nullable: false});
  const programInput = ui.input.choice('Program', {
    items: programs.map((p) => p.code), value: programs[0].code, nullable: false});
  const studyInput = ui.input.choice<string | null>('Study', {items: [null], value: null});
  const workstreamInput = ui.input.choice('Workstream', {
    items: [...WORKSTREAMS], value: 'modeling'});
  let linkedCodes: string[] = [];
  const compoundsInput = compoundsChipInput(() => linkedCodes);
  const tagsInput = tagsChipInput();
  const pathInput = ui.input.string('Path', {value: '', nullable: true});
  pathInput.setTooltip('Folder-like grouping, e.g. PK/exposure');
  const descriptionInput = ui.input.textArea('Description', {value: ''});
  const hint = ui.divText('', 'aa-republish-hint');

  const currentProgram = () => programs.find((p) => p.code === programInput.value)!;

  const refreshStudies = async () => {
    const studies = await grok.dapi.domains.table(T_STUDY).query({
      filter: [{property: 'program_id', operator: '=', value: currentProgram().id}], sort: 'protocol_code'});
    (studyInput as any).items = [null, ...studies.map((s) => s.protocol_code)];
    studyInput.value = null;
  };

  const refreshProgramCompounds = async () => {
    linkedCodes = await programCompoundCodes(currentProgram().id);
  };

  const refreshHint = async () => {
    const program = currentProgram();
    const study = studyInput.value;
    const name = nameInput.value?.trim();
    if (!name) {
      hint.textContent = '';
      return;
    }
    const live = await grok.dapi.domains.table(T_ALIGNMENT).query({filter: [
      {property: 'program_id', operator: '=', value: program.id},
      {property: 'study_id', operator: study == null ? 'is' : '=',
        value: study == null ? null : (await studyIdOf(program.id, study))},
      {property: 'name', operator: '=', value: name},
    ]});
    if (live.length > 0) {
      const next = Math.max(...live.map((r: any) => r.revision)) + 1;
      hint.textContent = `A publication named "${name}" exists in ${program.code}` +
        `${study ? ` / ${study}` : ''} — this publishes version ${next}; ` +
        'the current version stays visible until this one is approved.';
    } else
      hint.textContent = `This creates a new publication in ${program.code}.`;
  };

  const studyIdOf = async (programId: string, protocolCode: string | null): Promise<string | null> => {
    if (protocolCode == null)
      return null;
    const study = await grok.dapi.domains.table(T_STUDY).getByKey({protocol_code: protocolCode});
    return study?.id ?? null;
  };

  programInput.onChanged.subscribe(() => {
    refreshStudies().then(refreshHint);
    void refreshProgramCompounds();
  });
  studyInput.onChanged.subscribe(() => refreshHint());
  nameInput.onChanged.subscribe(() => refreshHint());
  await Promise.all([refreshStudies(), refreshProgramCompounds()]);
  await refreshHint();

  let dialog = ui.dialog('Publish to program');
  for (const input of [nameInput, programInput, studyInput, workstreamInput,
    compoundsInput, tagsInput, pathInput, descriptionInput])
    dialog = dialog.add(input);
  dialog = dialog.add(hint)
    .onOK(async () => {
      const program = currentProgram();
      const progress = DG.TaskBarProgressIndicator.create('Publishing…');
      try {
        const result = await publishWorkflowRun({
          ...(typeof source === 'string' ? {sourceMetaCallId: source} : {sourceCall: source}),
          programId: program.id,
          studyId: await studyIdOf(program.id, studyInput.value),
          name: nameInput.value!.trim(),
          workstream: workstreamInput.value ?? undefined,
          compounds: chips(compoundsInput.value),
          tags: chips(tagsInput.value),
          path: pathInput.value?.trim() || undefined,
          description: descriptionInput.value?.trim() || undefined,
        });
        grok.shell.info(`Published "${nameInput.value}" v${result.revision} (${result.status})`);
      } catch (e: any) {
        grok.shell.error(`Publish failed: ${e?.message ?? e}`);
      } finally {
        progress.close();
      }
    })
    .show({center: true, width: 500});

  const okButton = dialog.getButton('OK');
  const refreshOk = () => okButton.disabled = !nameInput.value?.trim();
  refreshOk();
  nameInput.onChanged.subscribe(refreshOk);
  return dialog;
}
