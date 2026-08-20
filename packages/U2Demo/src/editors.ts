import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {dartInputFor, inputForProperty, PropertyLike} from '@datagrok-libraries/u2/src/dg/index.js';

/** Property types a u2 editor can take over. The platform matches a `valueEditor` by
 * `options.propertyType`, and for these five that string is also the `DG.TYPE` the bridge has to
 * report — `dartInputFor` builds nothing until the platform asks for the element, so it has no
 * input to infer the type from. */
export const EDITOR_TYPES: string[] = [DG.TYPE.INT, DG.TYPE.FLOAT, DG.TYPE.BOOL,
  DG.TYPE.DATE_TIME, DG.TYPE.STRING];

const STORAGE_KEY = 'u2.valueEditors';
const registered = new Set<string>();

/** The types enabled in this browser profile. Nothing is enabled until someone opts in: a
 * registration takes over every matching property in the session, and only a reload undoes it. */
export function enabledEditorTypes(): string[] {
  return (localStorage.getItem(STORAGE_KEY) ?? '').split(',')
    .map((t) => t.trim()).filter((t) => EDITOR_TYPES.includes(t));
}

/** The types whose editors this session actually registered — a flag ticked since the last reload
 * is enabled but not yet active. */
export function registeredEditorTypes(): string[] {
  return [...registered];
}

export function setEditorTypes(types: string[]): void {
  localStorage.setItem(STORAGE_KEY, EDITOR_TYPES.filter((t) => types.includes(t)).join(','));
}

/** Registers one `role: valueEditor` func per type, so unmodified Dart forms, func-param dialogs
 * and property panels resolve u2 editors for properties of that type — the seam is consulted at
 * `input_base.dart:675`, before every type default and after an explicit `inputType`. The
 * platform replaces a func on re-registration; the guard keeps it to one per session. Rollback is
 * a reload: nothing on the JS side unregisters a func. */
export function registerU2ValueEditors(types: string[]): void {
  for (const type of types) {
    if (registered.has(type) || !EDITOR_TYPES.includes(type))
      continue;
    registered.add(type);
    grok.functions.register({
      signature: `object u2${type[0].toUpperCase()}${type.slice(1)}ValueEditor()`,
      tags: 'valueEditor',
      options: {propertyType: type},
      // assumeWritable: a func parameter carries no setter (`FuncParam` never assigns
      // `Property.set`), and the value belongs to the platform proxy rather than to the property
      run: () => dartInputFor((prop) => inputForProperty(prop as PropertyLike | null,
        {assumeWritable: true}), {dataType: type}),
    });
  }
}

export function registerEnabledEditors(): void {
  registerU2ValueEditors(enabledEditorTypes());
}
