import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {T_ALIGNMENT, T_PROGRAM, T_STUDY} from '../domain/constants';
import {publishWorkflowRun} from '../service/publication-service';

interface ProgramChoice {
  id: string;
  code: string;
  name: string;
}

/** The "Publish to program" dialog over a saved Compute2 workflow run (its meta
 * FuncCall id). The republish key is (program, study, name): a key match publishes
 * the next version of that publication, a new key creates a new publication. */
export async function showPublishDialog(sourceMetaCallId: string, defaultName?: string): Promise<void> {
  const programs: ProgramChoice[] =
    await grok.dapi.domains.table(T_PROGRAM).query({sort: 'code', limit: 1000});
  if (programs.length === 0) {
    grok.shell.warning('No programs are visible to you — ask an administrator to provision one');
    return;
  }

  const nameInput = ui.input.string('Name', {value: defaultName ?? ''});
  const programInput = ui.input.choice('Program', {
    items: programs.map((p) => p.code), value: programs[0].code, nullable: false});
  const studyInput = ui.input.choice<string | null>('Study', {items: [null], value: null});
  const workstreamInput = ui.input.choice('Workstream', {
    items: ['clinical', 'modeling', 'discovery', 'cmc', 'other'], value: 'modeling'});
  const compoundsInput = ui.input.string('Molecules', {value: '', nullable: true});
  compoundsInput.setTooltip('Comma-separated registration codes');
  const tagsInput = ui.input.string('Tags', {value: '', nullable: true});
  tagsInput.setTooltip('Comma-separated tags');
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
  });
  studyInput.onChanged.subscribe(() => refreshHint());
  nameInput.onChanged.subscribe(() => refreshHint());
  await refreshStudies();
  await refreshHint();

  const parseList = (value: string | null) =>
    (value ?? '').split(',').map((s) => s.trim()).filter((s) => s.length > 0);

  ui.dialog('Publish to program')
    .add(ui.divV([
      nameInput, programInput, studyInput, workstreamInput,
      compoundsInput, tagsInput, pathInput, descriptionInput, hint,
    ]))
    .onOK(async () => {
      const program = currentProgram();
      const progress = DG.TaskBarProgressIndicator.create('Publishing…');
      try {
        const result = await publishWorkflowRun({
          sourceMetaCallId,
          programId: program.id,
          studyId: await studyIdOf(program.id, studyInput.value),
          name: nameInput.value!.trim(),
          workstream: workstreamInput.value ?? undefined,
          compounds: parseList(compoundsInput.value),
          tags: parseList(tagsInput.value),
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
}
