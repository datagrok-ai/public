import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {historyUtils} from '@datagrok-libraries/compute-utils';
import {saveInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {T_ALIGNMENT, T_HISTORY, T_PROGRAM, T_STUDY} from '../domain/constants';
import {ensureProgram, ProgramInfo} from '../domain/security';

export const STEP_NQ = 'ArtifactAlignment:AaTestStep';
export const PROVIDER_NQ = 'ArtifactAlignment:AaTestProvider';

export interface SavedRunFixture {
  metaCallId: string;
  stepCallId: string;
}

/** A genuine saved compute2-style workflow run, built headlessly: one executed step
 * FuncCall with a dataframe output plus a meta call carrying the serialized pipeline
 * state that references it. */
export async function makeSavedRun(a: number = 3): Promise<SavedRunFixture> {
  const stepFunc = DG.Func.byName(STEP_NQ);
  const fc = stepFunc.prepare({a});
  await fc.call(false, undefined, {processed: true, report: false});
  fc.newId();
  const savedStep = await historyUtils.saveRun(fc);
  const state = {
    type: 'static',
    configId: 'root',
    uuid: crypto.randomUUID(),
    isReadonly: false,
    nqName: PROVIDER_NQ,
    friendlyName: 'AA fixture workflow',
    steps: [{
      type: 'funccall',
      configId: 'step1',
      uuid: crypto.randomUUID(),
      isReadonly: false,
      nqName: STEP_NQ,
      funcCallId: savedStep.id,
      friendlyName: 'Step 1',
    }],
  };
  const metaCall = await saveInstanceState(PROVIDER_NQ, state, {title: 'AA fixture run'});
  return {metaCallId: metaCall.id, stepCallId: savedStep.id};
}

const createdPrograms: ProgramInfo[] = [];

export async function makeTestProgram(): Promise<ProgramInfo> {
  const code = `ZZT-${Date.now()}-${Math.floor(Math.random() * 1e5)}`;
  const info = await ensureProgram({code, name: `Test program ${code}`});
  createdPrograms.push(info);
  return info;
}

export async function makeTestStudy(programId: string): Promise<{id: string, protocolCode: string}> {
  const protocolCode = `ZZS-${Date.now()}-${Math.floor(Math.random() * 1e5)}`;
  const [ins] = await grok.dapi.domains.table(T_STUDY).insert(
    {program_id: programId, protocol_code: protocolCode, phase: 'ph2', status: 'ongoing'});
  return {id: ins.id, protocolCode};
}

/** Deletes every run of the fixture functions together with the dataframes they
 * uploaded. Publishing freezes copies of these same functions' runs, so this sweeps
 * the frozen artifacts too. Time-boxed so a pre-existing backlog cannot blow the
 * category timeout — leftovers drain on subsequent runs. */
export async function cleanupTestRuns(budgetMs: number = 45000): Promise<void> {
  const started = Date.now();
  for (const nqName of [PROVIDER_NQ, STEP_NQ]) {
    const funcName = nqName.split(':')[1];
    while (Date.now() - started < budgetMs) {
      const page = await grok.dapi.functions.calls.allPackageVersions()
        .filter(`func.name="${funcName}"`).include('inputs, outputs').list({pageSize: 50});
      if (page.length === 0)
        break;
      await Promise.all(page.map(async (fc) => {
        // non-materialized dataframe params hold the uploaded table id
        const tableIds: string[] = [];
        for (const p of fc.inputParams.values() as Iterable<DG.FuncCallParam>) {
          if (p.property.propertyType === DG.TYPE.DATA_FRAME && typeof fc.inputs[p.name] === 'string')
            tableIds.push(fc.inputs[p.name]);
        }
        for (const p of fc.outputParams.values() as Iterable<DG.FuncCallParam>) {
          if (p.property.propertyType === DG.TYPE.DATA_FRAME && typeof fc.outputs[p.name] === 'string')
            tableIds.push(fc.outputs[p.name]);
        }
        for (const id of tableIds) {
          try {
            const table = await grok.dapi.tables.find(id);
            if (table != null)
              await grok.dapi.tables.delete(table);
          } catch (_) {/* already gone */}
        }
        // the plain calls source silently no-ops the delete — allPackageVersions is required
        await grok.dapi.functions.calls.allPackageVersions().delete(fc);
      }));
    }
  }
}

/** Removes everything the fixtures created: runs and their dataframes, version rows,
 * history, studies, program rows, and the per-program groups. */
export async function cleanupTestPrograms(): Promise<void> {
  await cleanupTestRuns();
  for (const program of createdPrograms.splice(0)) {
    const byProgram = DG.cond('program_id', '=', program.id);
    while ((await grok.dapi.domains.table(T_ALIGNMENT).deleteWhere(byProgram)).hasMore);
    while ((await grok.dapi.domains.table(T_HISTORY).deleteWhere(byProgram)).hasMore);
    while ((await grok.dapi.domains.table(T_STUDY).deleteWhere(byProgram)).hasMore);
    await grok.dapi.domains.table(T_PROGRAM).delete(program.id);
    for (const group of [program.viewers, program.contributors, program.approvers]) {
      try {
        await grok.dapi.groups.delete(group);
      } catch (_) {/* platform may refuse deleting groups still referenced */}
    }
  }
}
