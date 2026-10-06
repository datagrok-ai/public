/* Populate-from-a-saved-run: values are COPIED into the form's current call, so the call's
   identity (and absent outputs) is preserved — the Dart applyHistory semantics, server-backed.
   A stored run's dataframe params come back as TableInfo id strings and FILE params as
   `{id, name}` JSON (the compute-utils history convention) — both are materialized before the
   write. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import type {FuncCallForm} from './func-form.js';

export async function materializeInputs(picked: DG.FuncCall): Promise<Map<string, any>> {
  const out = new Map<string, any>();
  for (const p of picked.inputParams.values() as Iterable<DG.FuncCallParam>) {
    let v = picked.inputs[p.name];
    try {
      if (typeof v === 'string') {
        if (p.property.propertyType === DG.TYPE.DATA_FRAME)
          v = await grok.dapi.tables.getTable(v);
        else if (p.property.propertyType === DG.TYPE.FILE) {
          const {id, name} = JSON.parse(v);
          const fileInfo = DG.FileInfo.fromBytes(name, await grok.dapi.files.readAsBytes(id));
          fileInfo.id = id;
          v = fileInfo;
        }
      }
      out.set(p.name, v);
    } catch (_) { /* the Dart _materialize tolerance: a vanished table skips the one param */ }
  }
  return out;
}

/** Copies the picked run's input values into the form's CURRENT call. Only params the current
 * call carries are written — a stale pick from another function is inherently a no-op. */
export async function applyHistory(form: FuncCallForm, picked: DG.FuncCall): Promise<void> {
  const values = await materializeInputs(picked);
  const current = new Set([...form.source.inputParams.values()].map((p) => p.name));
  for (const [name, v] of values) {
    if (current.has(name) && v !== undefined)
      form.source.setParamValue(name, v); // two-way binding refreshes the inputs (echo-suppressed)
  }
}

export async function applyHistoryById(form: FuncCallForm, id: string): Promise<void> {
  const picked = await grok.dapi.functions.calls.allPackageVersions()
    .include('func.package,func.params,inputs').find(id);
  if (picked != null)
    await applyHistory(form, picked);
}
