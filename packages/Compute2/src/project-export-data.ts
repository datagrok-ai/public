import * as DG from 'datagrok-api/dg';

export interface DfExportEntry {
  name: string;
  isInput: boolean;
  df: DG.DataFrame;
  viewers: {type: string, options: Record<string, any>}[];
}

const NON_SYNCABLE_INPUT_TYPES = new Set<string>([
  DG.TYPE.DATA_FRAME_LIST, DG.TYPE.FILE, DG.TYPE.BLOB,
  DG.TYPE.COLUMN, DG.TYPE.COLUMN_LIST, DG.TYPE.GRAPHICS, DG.TYPE.OBJECT, DG.TYPE.DYNAMIC,
]);

/** The producing call serialized by the platform itself (`prepare().toString()`), or null
 *  when the call cannot be replayed — outputs then publish as static snapshots.
 *  A dataframe input serializes as its quoted table name (same as the core's derived-table
 *  provenance, resolved among project children on open), so the exported clone — carrying
 *  the final deduped name — is what goes into prepare(); a df input that is not exported
 *  would leave a dangling reference, so it disables sync. Unset inputs replay as nulls and
 *  do not affect syncability. */
function serializeProducingCall(call: DG.FuncCall, entries: DfExportEntry[]): string | null {
  const valued = [...call.inputParams.values()]
    .map((p) => p.property)
    .filter((p) => call.inputs[p.name] != null);
  if (valued.some((p) => NON_SYNCABLE_INPUT_TYPES.has(p.propertyType as string)))
    return null;
  try {
    const inputs: Record<string, any> = {};
    for (const prop of valued) {
      let value = call.inputs[prop.name];
      if (prop.propertyType === DG.TYPE.DATA_FRAME) {
        const entry = entries.find((e) => e.isInput && e.name === prop.name);
        if (!entry)
          return null;
        value = entry.df;
      }
      inputs[prop.name] = value;
    }
    return call.func.prepare(inputs).toString();
  } catch (e) {
    console.warn('Compute2: could not serialize the producing call', e);
    return null;
  }
}

/** Stamp each output clone with the producing call (`.script` + `.VariableName` tags) so the
 *  core Save-project dialog offers data sync — identical producing calls dedup on project open,
 *  so the function runs once and every table binds its own output. */
function stampCreationScripts(call: DG.FuncCall, entries: DfExportEntry[]) {
  const outputs = entries.filter((e) => !e.isInput);
  if (outputs.length === 0)
    return;
  const callStr = serializeProducingCall(call, entries);
  if (callStr == null)
    return;
  const accessor = [...call.outputParams.values()].length > 1;
  let ts = Date.now();
  for (const entry of outputs) {
    entry.df.setTag(DG.Tags.VariableName, entry.name);
    entry.df.setTag(DG.Tags.CreationScript,
      `${entry.name} = ${callStr}${accessor ? '.' + entry.name : ''} //{"timestamp": ${ts++}}`);
  }
}

/** Clones each entry's dataframe under a unique non-empty name (a null or blank df name
 *  falls back to the parameter name — an unnamed df would otherwise serialize as a `null`
 *  argument in the creation script) and stamps output clones for data sync. */
export function prepareExportEntries(call: DG.FuncCall, entries: DfExportEntry[]): DfExportEntry[] {
  const usedNames = new Set<string>();
  const cloned = entries.map((entry) => {
    const df = entry.df.clone();
    const baseName = (entry.df.name ?? '').trim() || entry.name;
    let name = baseName;
    for (let i = 2; usedNames.has(name); i++)
      name = `${baseName} (${i})`;
    usedNames.add(name);
    df.name = name;
    return {...entry, df};
  });
  stampCreationScripts(call, cloned);
  return cloned;
}
