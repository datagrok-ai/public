import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {IControllerBase} from '../RuntimeControllers';
import {RuleSource} from '../config/PipelineConfiguration';

export type ValidatorVerdict = {message: string, isError: boolean, isHelper: boolean};

let warnedNoEval = false;

// the annotation's own validators, resolved and run by the platform (1.28+)
function annotationVerdicts(
  controller: IControllerBase, spec: {input: string, call?: string},
): ValidatorVerdict[] | Promise<ValidatorVerdict[]> {
  if (!spec.call || !controller.hasCall(spec.call))
    return [];
  const call = controller.getFirst<DG.FuncCall | undefined>(spec.call);
  if (!call)
    return [];
  if (typeof (call as any).evalParamValidators !== 'function') {
    if (!warnedNoEval) {
      warnedNoEval = true;
      console.warn('RTD: FuncCall.evalParamValidators is not available on this platform, named validators are skipped');
    }
    return [];
  }
  const io = controller.getMatchedPositions(spec.input)[0]?.ioName;
  if (!io)
    return [];
  return call.evalParamValidators(io).catch((e) => [failed('', e)]);
}

function failed(name: string, e: unknown): ValidatorVerdict {
  const detail = e instanceof Error ? e.message : String(e);
  return {message: `Couldn't validate (${name ? name + ': ' : ''}${detail})`, isError: false, isHelper: false};
}

// the documented validator contract: one argument, null/true valid, string is the
// message, false is a generic failure naming the function
async function callValidator(name: string, value: any): Promise<ValidatorVerdict | undefined> {
  const func = DG.Func.byName(name);
  if (!func)
    return failed(name, 'function not found');
  const param = func.inputs[0]?.name;
  if (!param)
    return failed(name, 'validator takes no input');
  try {
    const result = await grok.functions.call(name, {[param]: value});
    if (result == null || result === true)
      return undefined;
    if (result === false)
      return {message: `Validation failed: ${name}`, isError: true, isHelper: false};
    return {message: String(result), isError: true, isHelper: false};
  } catch (e) {
    return failed(name, e);
  }
}

function namedVerdicts(value: any, names: string[]): Promise<ValidatorVerdict[]> {
  return Promise.all(names.map((name) => callValidator(name, value)))
    .then((verdicts) => verdicts.filter((verdict): verdict is ValidatorVerdict => !!verdict));
}

// synchronous whenever nothing is called (absent value, mock mode, old platform), so those
// links stay batchable and virtual-time tests see their results immediately
function resolveValidators(
  controller: IControllerBase, spec: RuleSource['validators'],
): ValidatorVerdict[] | Promise<ValidatorVerdict[]> {
  if (spec.names) {
    const value = controller.getFirst(spec.input);
    if (value == null || value === '')
      return [];
    return namedVerdicts(value, spec.names);
  }
  return annotationVerdicts(controller, spec);
}

/** Resolves the values a rule declares in `sources`; each alias becomes a context variable.
 *  Returns a plain object when nothing had to be awaited. */
export function resolveSources(
  controller: IControllerBase, sources: Record<string, RuleSource> | undefined,
): Record<string, any> | Promise<Record<string, any>> {
  const resolved: Record<string, any> = {};
  const pending: Promise<void>[] = [];
  for (const [alias, source] of Object.entries(sources ?? {})) {
    if (!('validators' in source))
      throw new Error(`Unknown rule source ${JSON.stringify(source)} for alias ${alias}`);
    const value = resolveValidators(controller, source.validators);
    if (value instanceof Promise)
      pending.push(value.then((verdicts) => { resolved[alias] = verdicts; }));
    else
      resolved[alias] = value;
  }
  return pending.length ? Promise.all(pending).then(() => resolved) : resolved;
}
