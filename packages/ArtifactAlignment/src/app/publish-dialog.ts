import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {T_ALIGNMENT, T_COMPOUND, T_PROGRAM, T_PROGRAM_COMPOUND, T_STUDY, T_TAG} from '../domain/constants';
import {publishWorkflowRun} from '../service/publication-service';

interface ProgramChoice {
  id: string;
  code: string;
  name: string;
}

/** The "Publish to program" dialog over a Compute2 run — a saved run's meta FuncCall
 * id, or a live in-memory FuncCall (published without a prior history save). The
 * republish key is (program, study, name): a key match publishes the next version
 * of that publication, a new key creates a new publication. */
export async function showPublishDialog(source: string | DG.FuncCall, defaultName?: string): Promise<void> {
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
  // chips with registry-backed suggestions; empty text suggests the program's own
  // compounds. allowNew is off — a typo must not mint a new registry compound.
  let programCompoundCodes: string[] = [];
  const compoundsInput = ui.input.tags('Molecules', {
    allowNew: false,
    getSuggestions: async (text: string) => {
      const t = text?.trim() ?? '';
      if (t === '')
        return programCompoundCodes;
      const rows = await grok.dapi.domains.table(T_COMPOUND).query({
        filter: DG.cond('registration_code', 'like', `%${t}%`), sort: 'registration_code', limit: 20});
      return rows.map((r: any) => r.registration_code);
    },
  });
  compoundsInput.setTooltip('Registration codes from the compound registry');
  const tagsInput = ui.input.tags('Tags', {
    allowNew: true,
    createNewItem: (text: string) => text.trim(),
    getSuggestions: async (text: string) => {
      const rows = await grok.dapi.domains.table(T_TAG).query({
        filter: DG.cond('name', 'like', `%${text?.trim() ?? ''}%`), sort: 'name', limit: 20});
      return rows.map((r: any) => r.name);
    },
  });
  tagsInput.setTooltip('Existing tags are suggested; a new name creates the tag');
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
    const links = await grok.dapi.domains.table(T_PROGRAM_COMPOUND).query({
      filter: DG.cond('program_id', '=', currentProgram().id), limit: 100});
    programCompoundCodes = links.length === 0 ? [] :
      (await grok.dapi.domains.table(T_COMPOUND).query({
        filter: DG.or(...links.map((l: any) => DG.cond('id', '=', l.compound_id))), limit: 100}))
        .map((c: any) => c.registration_code).sort();
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

  const chips = (value: unknown): string[] =>
    ((value ?? []) as string[]).map((s) => `${s}`.trim()).filter((s) => s.length > 0);

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
}
