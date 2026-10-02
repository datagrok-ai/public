export type {ObservableLike, IEventType, IRectBounds, IInputStatus, IWidgetStatus, FuncLike,
  NamedProperty, BindProp, BindSource, ComponentMetaBase} from 'datagrok-api/u2core';
import type {IInputStatus} from 'datagrok-api/u2core';
import type {FieldAccess} from './access.js';

/** An input's status with the access level a schema-driven form applied to it — why a field is
 * text (`readonly`) or an editor (`editable`); a hidden field has no status at all. */
export interface IFieldStatus extends IInputStatus {
  access: FieldAccess;
}

/** The {@link BindProp.type} label inferred from a live value — what the binding picker shows. */
export function bindTypeOf(value: unknown): string {
  if (typeof value === 'string')
    return 'string';
  if (typeof value === 'number')
    return Number.isInteger(value) ? 'int' : 'double';
  if (typeof value === 'boolean')
    return 'bool';
  if (Array.isArray(value) && value.every((item) => typeof item === 'string'))
    return 'string_list';
  return 'object';
}
