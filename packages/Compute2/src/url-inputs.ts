import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

export interface UrlInputsResult {
  patch: Map<string, any>;
  entityIds: Map<string, string>;
  warnings: string[];
}

// 'id' is the run-loading param handled by RFVApp; platform startup params (q, layout, cmd, ...)
// are not reserved — the platform acts on them but leaves them in startUri, and an input
// sharing such a name should still be settable
const RESERVED_PARAMS = new Set(['id']);

const SCALAR_TYPES = new Set<string>([DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BOOL, DG.TYPE.STRING, DG.TYPE.DATE_TIME]);

function parseScalar(type: string, raw: string): {ok: boolean, value?: any} {
  if (raw.trim() === '' && type !== DG.TYPE.STRING)
    return {ok: false};
  switch (type) {
  case DG.TYPE.INT: {
    const num = Number(raw);
    return Number.isInteger(num) ? {ok: true, value: num} : {ok: false};
  }
  case DG.TYPE.FLOAT: {
    const num = Number(raw);
    return Number.isFinite(num) ? {ok: true, value: num} : {ok: false};
  }
  case DG.TYPE.BOOL: {
    const val = raw.toLowerCase();
    if (val === 'true' || val === '1') return {ok: true, value: true};
    if (val === 'false' || val === '0') return {ok: true, value: false};
    return {ok: false};
  }
  case DG.TYPE.STRING:
    return {ok: true, value: raw};
  case DG.TYPE.DATE_TIME: {
    const date = dayjs(raw);
    return date.isValid() ? {ok: true, value: date} : {ok: false};
  }
  default:
    return {ok: false};
  }
}

/** Matches URL search params to the call's inputs and pre-loads entity-valued ones
 *  (dataframe/file params carry entity ids), so the patch can later be applied synchronously. */
export async function parseUrlInputs(call: DG.FuncCall, urlParams: URLSearchParams): Promise<UrlInputsResult> {
  const patch = new Map<string, any>();
  const entityIds = new Map<string, string>();
  const warnings: string[] = [];
  const byName = new Map([...call.inputParams.values()].map((p) => [p.property.name, p]));

  for (const [name, raw] of urlParams.entries()) {
    if (RESERVED_PARAMS.has(name))
      continue;
    // unmatched params are skipped silently: platform params (q, layout, browse, ...) stay in the URL
    const param = byName.get(name);
    if (!param)
      continue;
    const type = param.property.propertyType as string;
    if (type === DG.TYPE.DATA_FRAME) {
      try {
        patch.set(name, await grok.dapi.tables.getTable(raw));
        entityIds.set(name, raw);
      } catch {
        warnings.push(`Failed to load table "${raw}" for input "${name}"`);
      }
    } else if (type === DG.TYPE.FILE) {
      try {
        const entity = await grok.dapi.entities.find(raw);
        if (entity instanceof DG.FileInfo) {
          patch.set(name, entity);
          entityIds.set(name, raw);
        } else
          warnings.push(`Entity "${raw}" for input "${name}" is not a file`);
      } catch {
        warnings.push(`Failed to load file "${raw}" for input "${name}"`);
      }
    } else if (SCALAR_TYPES.has(type)) {
      const res = parseScalar(type, raw);
      if (res.ok)
        patch.set(name, res.value);
      else
        warnings.push(`Cannot parse "${raw}" as ${type} for input "${name}"`);
    } else
      warnings.push(`Input "${name}" of type ${type} is not supported in URL`);
  }
  return {patch, entityIds, warnings};
}

export function applyUrlInputs(call: DG.FuncCall, patch: Map<string, any>) {
  for (const [name, value] of patch)
    call.inputs[name] = value;
}

// `nullable` is not consulted: it defaults to true for every non-string type, so it carries
// no authored signal — `optional: true` is the marker for inputs that may stay empty
export function missingMandatoryInputs(call: DG.FuncCall): string[] {
  return [...call.inputParams.values()]
    .filter((p) => !p.property.isOptional && call.inputs[p.property.name] == null)
    .map((p) => p.property.name);
}

/** Current view URL with input values as search params; dataframe/file inputs are included
 *  only when their entity id is known (loaded from a URL), others are reported as skipped. */
export function buildInputsUrl(call: DG.FuncCall, entityIds: Map<string, string>): {url: string, skipped: string[]} {
  const url = new URL(window.location.href);
  url.search = '';
  url.hash = '';
  const skipped: string[] = [];
  for (const param of call.inputParams.values()) {
    const name = param.property.name;
    const value = call.inputs[name];
    if (value == null)
      continue;
    const type = param.property.propertyType as string;
    if (type === DG.TYPE.DATA_FRAME || type === DG.TYPE.FILE) {
      const id = entityIds.get(name);
      if (id)
        url.searchParams.set(name, id);
      else
        skipped.push(name);
    } else if (type === DG.TYPE.DATE_TIME)
      url.searchParams.set(name, dayjs(value).toISOString());
    else if (SCALAR_TYPES.has(type))
      url.searchParams.set(name, String(value));
    else
      skipped.push(name);
  }
  return {url: url.toString(), skipped};
}

export async function copyText(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    const ta = document.createElement('textarea');
    ta.value = text;
    ta.style.position = 'fixed';
    document.body.appendChild(ta);
    ta.select();
    try {
      return document.execCommand('copy');
    } catch {
      return false;
    } finally {
      ta.remove();
    }
  }
}
