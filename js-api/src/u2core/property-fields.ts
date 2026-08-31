// Type-only: u2core stays a value-level leaf — never add a runtime import here.
import type {IProperty} from '../entities/property-meta.js';

export type {IProperty, IPropertyMeta} from '../entities/property-meta.js';

const PROPERTY_FIELDS: (keyof IProperty)[] = ['name', 'caption', 'friendlyName', 'propertyType', 'type',
  'semType', 'description', 'category', 'choices', 'inputType', 'editor', 'nullable', 'min', 'max', 'step',
  'format', 'units', 'showSlider', 'showPlusMinus', 'get', 'set'];

/** A plain copy of the named fields — a spread of a `DG.Property` copies nothing but its `dart`
 * handle, since every field is a prototype getter. */
export function propertyFields(p: IProperty): IProperty {
  const out: Record<string, unknown> = {};
  for (const key of PROPERTY_FIELDS) {
    if (p[key] !== undefined)
      out[key] = p[key];
  }
  return out as unknown as IProperty;
}
