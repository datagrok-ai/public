/**
 * Canonical property metadata interfaces.
 * u2core type-imports this file — imports here must stay `import type` at the value level;
 * never add a runtime import.
 * @module entities/property-meta
 */

import type {PropertyGetter, PropertySetter, ValueValidator} from './types';
import type {ColumnTypeFilter} from '../const';

/** The pure property-metadata vocabulary: what a form generator, a property grid or a spec
 * registry reads off a property descriptor. See also {@link IProperty}, {@link Property}. */
export interface IPropertyMeta {

  /** Property name */
  name?: string;

  /** Property data type. See {@link TYPE}. */
  type?: string;

  /** Property data type; same as {@link type}, and read first where both appear (as on a `DG.Property`). */
  propertyType?: string;

  /** Platform `InputType` ('Radio', 'Color', 'Tags'…): the editor named here wins over every
   * other hint, as it does in `input_base.dart:668`. */
  inputType?: string;

  /** Whether an empty value is allowed. This is used by validators. */
  nullable?: boolean;

  /** Property description */
  description?: string;

  /** Semantic type */
  semType?: string;

  /** Units of measurement. See also: [postfix] */
  units?: string;

  /** Minimum value. Applicable to numerical properties only */
  min?: number;

  /** Maximum value. Applicable to numerical properties only */
  max?: number;

  /** Step to be used in a slider. Only applies to numerical properties. */
  step?: number;

  /** Whether a slider appears next to the number input. Applies to numerical columns only. */
  showSlider?: boolean;

  /** Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only. */
  showPlusMinus?: boolean;

  /** List of choices. Applicable to string properties only */
  choices?: string[];

  /** Platform editor hint ('textarea', 'password', 'switch', 'slider'), consulted next, and only
   * for the types Dart honors it for (`input_base.dart:702-728`). */
  editor?: string;

  /** Corresponding category on the context panel; groups the property in a property panel. */
  category?: string;

  /** Value format, such as '0.000', applied through `DG.format` where the platform is loaded. */
  format?: string;

  /** Custom field friendly name shown in [PropertyGrid] */
  friendlyName?: string;

  /** Custom field caption shown in [PropertyGrid]
   * @deprecated The property will be removed soon. Use {@link friendlyName} instead */
  caption?: string;
}

/** Represents a property: the metadata plus accessors and behavior options.
 * See also {@link Property}. */
export interface IProperty extends IPropertyMeta {

  /** Reads the property value off an object. A property with no {@link set} renders read-only. */
  get?: PropertyGetter;

  /** Writes the property value to an object. */
  set?: PropertySetter;

  /** Initial value used when initializing UI. See also {@link defaultValue} */
  initialValue?: any;

  /** Default value used for deserialization and cloning. See also {@link initialValue}. */
  defaultValue?: any;

  /** List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names.
   * Signature: validator(x: DG.Type): string | null.
   * [null] indicates that the value is valid, [string] describes a validation error. */
  validators?: string[];

  /** List of value validators (functions that take a value and return error message or null) */
  valueValidators?: ValueValidator<any>[];

  /** Name of the corresponding JavaScript field. No need to specify it if it is the same as name. */
  fieldName?: string;

  tags?: any;

  /** Additional options. */
  options?: any;

  /** Filter for columns, can be numerical, numerical_no_datetime, categorical, datetime,
   * categorical_or_datetime, or directly a column type (string, int...)
   * Applicable when type = Column */
  columnTypeFilter?: ColumnTypeFilter | null;

  viewer?: string;

  /** Whether the property should be editable via the UI */
  userEditable?: boolean;
}
