/** Saving a run to a program from the app: the runs of catalogued models this session made, and
 * what the save-to-program dialog offers for one of them. The rules live in ./model/save-run. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb} from './generated/db';
import {readRunAnnotations} from './model/model-runs';
import {Choice, programChoices, SaveRunChoices} from './model/save-run';
import {splitCodes} from './model/catalogue';
import {portfolioByCode} from './portal';
import {filterValue, ROW_LIMIT} from './rows';

/** The codes of the catalogued models, read once per session; refreshed by `catalogueCodes(true)`. */
let codes: Promise<Set<string>> | null = null;
export function catalogueCodes(refresh = false): Promise<Set<string>> {
  if (!codes || refresh) {
    codes = portalDb.models.query({limit: ROW_LIMIT}).then((ms) => new Set(ms.map((m) => m.code)))
      .catch(() => {
        // not cached: one failed read must not turn run tracking off for the session
        codes = null;
        return new Set<string>();
      });
  }
  return codes;
}

const MAX_SESSION_RUNS = 20;
const runs: DG.FuncCall[] = [];

/** The runs of catalogued models made in this session since the package loaded, newest first.
 * A run made before the package loaded is not here: save it with `saveRunToProgram(call)`. */
export function sessionRuns(): DG.FuncCall[] {
  return [...runs];
}

let tracking = false;

/** Starts keeping the runs of catalogued models, and calls [onRun] for each. Once per session. */
export function trackModelRuns(onRun: (call: DG.FuncCall) => void): void {
  if (tracking)
    return;
  tracking = true;
  grok.functions.onAfterRunAction.subscribe(async (call: DG.FuncCall) => {
    const name = call?.func?.nqName;
    if (!name || !(await catalogueCodes()).has(name) || runs.includes(call))
      return;
    runs.unshift(call);
    runs.splice(MAX_SESSION_RUNS);
    onRun(call);
  });
}

/** The programs of a portfolio the user can see, as the dialog lists them. */
export async function programsOf(portfolioCode: string): Promise<Choice[]> {
  const pf = await portfolioByCode(portfolioCode);
  return programChoices(await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)}`, limit: ROW_LIMIT}));
}

/** Everything the dialog offers for [call]: portfolios, unit operations, run types, and the
 * model's own catalogue entries to fill in from. */
export async function saveRunChoices(call: DG.FuncCall | null,
  context?: {portfolio?: string | null, program?: string | null} | null): Promise<SaveRunChoices> {
  const [portfolios, ops, types] = await Promise.all([
    portalDb.portfolios.query({limit: ROW_LIMIT}), portalDb.unitOperations.query({limit: ROW_LIMIT}), portalDb.runTypes.query({limit: ROW_LIMIT})]);
  const model = call?.func ? await portalDb.models.getByKey({code: call.func.nqName}).catch(() => null) : null;
  const runType = model?.default_run_type_id ? types.find((t) => t.id === model.default_run_type_id) : undefined;
  const byName = <T extends {code: string, name: string}>(xs: T[]) => xs.map((x) => ({code: x.code, name: x.name}))
    .sort((a, b) => a.name.localeCompare(b.name));
  return {
    portfolios: byName(portfolios),
    modelPortfolios: splitCodes(model?.portfolios),
    unitOperations: ops.sort((a, b) => (a.ordinal ?? 0) - (b.ordinal ?? 0) || a.name.localeCompare(b.name))
      .map((x) => ({code: x.code, name: x.name})),
    modelUnitOperations: splitCodes(model?.unit_operations),
    runTypes: byName(types),
    defaultRunType: runType?.code ?? null,
    annotations: call ? readRunAnnotations(call.options as any) : null,
    context: context ?? null,
  };
}
