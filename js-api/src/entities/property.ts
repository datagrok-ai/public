/**
 * Property classes for entity metadata.
 * @module entities/property
 */

import {ColumnType, SemType, TYPE, Type} from "../const";
import {toDart, toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import {IDartApi} from "../api/grok_api.g";
import {InputType} from "../api/d4.api.g";
import {PropertyGetter, PropertySetter, ValueValidator} from "./types";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Represents a property.
 * See also {@link Property}. */
export interface IProperty {

  /** Property name */
  name?: string;

  /** Property data type. See {@link TYPE}. */
  type?: string;

  /** Property input type */
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

  /** Initial value used when initializing UI. See also {@link defaultValue} */
  initialValue?: any;

  /** Default value used for deserialization and cloning. See also {@link initialValue}. */
  defaultValue?: any;

  /** Custom editor (such as slider or text area) */
  editor?: string;

  /** Corresponding category on the context panel */
  category?: string;

  /** Value format, such as '0.000' */
  format?: string;

  /** Whether the property should be editable via the UI */
  userEditable?: boolean;

  /** List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names.
   * Signature: validator(x: DG.Type): string | null.
   * [null] indicates that the value is valid, [string] describes a validation error. */
  validators?: string[];

  /** List of value validators (functions that take a value and return error message or null) */
  valueValidators?: ValueValidator<any>[];

  /** Custom field caption shown in [PropertyGrid]
   * @deprecated The property will be removed soon. Use {@link friendlyName} instead */
  caption?: string;

  /** Custom field friendly name shown in [PropertyGrid] */
  friendlyName?: string;

  /** Name of the corresponding JavaScript field. No need to specify it if it is the same as name. */
  fieldName?: string;

  tags?: any;

  /** Additional options. */
  options?: any;

  /** Filter for columns, can be numerical, categorical or directly a column type (string, int...)
   * Applicable when type = Column */
  columnTypeFilter?: ColumnType | 'numerical' | 'categorical' | null;


  viewer?: string;
}


/** Properties of properties. See also {@link Property.propertyOptions}. */
interface IPropertyMeta {
  applicableTo?: string;
}


/**
 * Strongly-typed property associated with an object.
 * Used for reflection, serialization, UI generation, and other introspection-dependent tasks.
 *
 * Samples:
 */
export class Property implements IProperty {
  public readonly dart: any;
  public options: any;

  constructor(dart: any) {
    this.dart = dart;
    this.options = new MapProxy(api.grok_Property_Get_Options(this.dart));
  }

  /** Property getter is a function that accepts one parameter (item)
   * and returns the property value.
   *
   * @returns {PropertyGetter} */
  get get(): PropertyGetter { return api.grok_Property_Get_Get(this.dart); }
  set get(x: PropertyGetter) { api.grok_Property_Set_Get(this.dart, x); }

  /** Property setter */
  get set(): PropertySetter { return api.grok_Property_Get_Set(this.dart); }
  set set(x: PropertySetter) { api.grok_Property_Set_Set(this.dart, x); }

  /** Property name */
  get name(): string { return api.grok_Property_Get_Name(this.dart); }
  set name(s: string) { api.grok_Property_Set_Name(this.dart, s); }

  /** Custom field caption shown in the UI
   * @deprecated The property will be removed soon. Use {@link friendlyName} instead */
  get caption(): string { return api.grok_Property_Get_Caption(this.dart); }
  set caption(s: string) { api.grok_Property_Set_Caption(this.dart, s); }

  /** Custom field caption shown in the UI */
  get friendlyName(): string { return api.grok_Property_Get_Caption(this.dart); }
  set friendlyName(s: string) { api.grok_Property_Set_Caption(this.dart, s); }

  /** Property category */
  get category(): string { return api.grok_Property_Get_Category(this.dart); }
  set category(s: string) { api.grok_Property_Set_Category(this.dart, s); }

  /** Property type. Same as {@link propertyType} */
  get type(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set type(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Property type. Same as {@link type} */
  get propertyType(): TYPE { return api.grok_Property_Get_PropertyType(this.dart); }
  set propertyType(s: TYPE) { api.grok_Property_Set_PropertyType(this.dart, s); }

  /** Property subtype */
  get propertySubType(): TYPE { return api.grok_Property_Get_PropertySubType(this.dart); }

  /** Applies to viewers properties whether to include the property in the layout or not. */
  get includeInLayout(): boolean { return api.grok_Property_Get_IncludeInLayout(this.dart); }
  set includeInLayout(s: boolean) { api.grok_Property_Set_IncludeInLayout(this.dart, s); }

  /** Semantic type */
  get semType(): SemType | string { return api.grok_Property_Get_SemType(this.dart); }
  set semType(s: SemType | string) { api.grok_Property_Set_SemType(this.dart, s); }

  /** Input type. See also {@link InputType} */
  get inputType(): string { return api.grok_Property_Get_InputType(this.dart); }
  set inputType(s: string) { api.grok_Property_Set_InputType(this.dart, s); }

  /** Description */
  get description(): string { return api.grok_Property_Get_Description(this.dart); }
  set description(s: string) { api.grok_Property_Set_Description(this.dart, s); }

  /** Nullable */
  get nullable(): boolean { return api.grok_Property_Get_Nullable(this.dart); }
  set nullable(s: boolean) { api.grok_Property_Set_Nullable(this.dart, s); }

  /** Initial value used when initializing UI */
  get initialValue(): string { return toJs(api.grok_Property_Get_InitialValue(this.dart)); }
  set initialValue(s: string) { api.grok_Property_Set_InitialValue(this.dart, toDart(s)); }

  /** Default value */
  get defaultValue(): any { return toJs(api.grok_Property_Get_DefaultValue(this.dart)); }
  set defaultValue(s: any) { api.grok_Property_Set_DefaultValue(this.dart, toDart(s)); }

  /** Property editor */
  get editor(): string { return api.grok_Property_Get(this.dart, 'editor'); }
  set editor(s: string) { api.grok_Property_Set(this.dart, 'editor', s); }

  /** Units of measurement. */
  get units(): string { return api.grok_Property_Get_Units(this.dart); }
  set units(s: string) { api.grok_Property_Set_Units(this.dart, s); }

  /** Format to be used for displaying the value. */
  get format(): string { return api.grok_Property_Get_Format(this.dart); }
  set format(s: string) { api.grok_Property_Set_Format(this.dart, s); }

  /** Whether a user can edit this property from the UI. */
  get userEditable(): boolean { return api.grok_Property_Get_UserEditable(this.dart); }
  set userEditable(s: boolean) { api.grok_Property_Set_UserEditable(this.dart, s); }

  /** Minimum value. Used when constructing UI (sliders), validating, etc. */
  get min(): number { return api.grok_Property_Get_Min(this.dart); }
  set min(s: number) { api.grok_Property_Set_Min(this.dart, s); }

  /** Maximum value. Used when constructing UI (sliders), validating, etc. */
  get max(): number { return api.grok_Property_Get_Max(this.dart); }
  set max(s: number) { api.grok_Property_Set_Max(this.dart, s); }

  /** Step to be used in a slider. Only applies to numerical properties. */
  get step(): number { return api.grok_Property_Get_Step(this.dart); }
  set step(s: number) { api.grok_Property_Set_Step(this.dart, s); }

  /** Whether a slider appears next to the number input. Applies to numerical columns only. */
  get showSlider(): boolean { return api.grok_Property_Get_ShowSlider(this.dart); }
  set showSlider(s: boolean) { api.grok_Property_Set_ShowSlider(this.dart, s); }

  /** Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only. */
  get showPlusMinus(): boolean { return api.grok_Property_Get_ShowPlusMinus(this.dart); }
  set showPlusMinus(s: boolean) { api.grok_Property_Set_ShowPlusMinus(this.dart, s); }

  /** List of possible values of that property.
   *  PropertyGrid will use it to populate combo boxes.
   *  @returns {Array<string>} */
  get choices(): string[] { return api.grok_Property_Get_Choices(this.dart); }
  set choices(x: string[]) { api.grok_Property_Set_Choices(this.dart, x); }

  /** Validation conditions to be checked when editing the property.
   * See also {@link PropertyValidator}. */
  get validators(): string[] { return api.grok_Property_Get_Validators(this.dart); }
  set validators(x: string[]) { api.grok_Property_Set_Validators(this.dart, x); }

  get isVectorizable(): boolean { return api.grok_Property_Get_IsVectorizable(this.dart); }

  get vectorName(): string { return api.grok_Property_Get_VectorName(this.dart); }

  /** Column type filter (previously "columnFilter") */
  get columnTypeFilter(): ColumnType | 'numerical' | 'categorical' | null {
    return api.grok_Property_Get_ColumnTypeFilter(this.dart);
  }

  /** Applies the specified options */
  fromOptions(opt?: IProperty): Property {
    if (opt)
      api.grok_Property_Options(this.dart, opt);
    return this;
  }

  /** Creates a property */
  static create(name: string, type: Type,
                getter: PropertyGetter,
                setter: PropertySetter,
                defaultValue: any = null): Property {
    return new Property(api.grok_Property(name, type, getter, setter, toDart(defaultValue)));
  }

  /** Creates an integer property */
  static int(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.INT, getter, setter, defaultValue);
  }

  /** Creates a float property */
  static float(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.FLOAT, getter, setter, defaultValue);
  }

  /** Creates a string property */
  static string(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.STRING, getter, setter, defaultValue);
  }

  /** Creates a bool property */
  static bool(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.BOOL, getter, setter, defaultValue);
  }

  /** Creates a datetime property */
  static dateTime(name: string, getter: PropertyGetter, setter: PropertySetter, defaultValue: any): Property {
    return Property.create(name, TYPE.DATE_TIME, getter, setter, defaultValue);
  }

  /** Creates property for the JavaScript objects with the corresponding property name */
  static js(name: string, type: TYPE, options?: IProperty): Property {
    return Property.create(name, type,
      (x: any) => x[name],
      function (x: any, v: any) { x[name] = v; },
      options?.defaultValue).fromOptions(options);
  }

  static jsInt(name: string, options?: IProperty): Property { return Property.js(name, TYPE.INT, options); }
  static jsBool(name: string, options?: IProperty): Property { return Property.js(name, TYPE.BOOL, options); }
  static jsFloat(name: string, options?: IProperty): Property { return Property.js(name, TYPE.FLOAT, options); }
  static jsString(name: string, options?: IProperty): Property { return Property.js(name, TYPE.STRING, options); }
  static jsDateTime(name: string, options?: IProperty): Property { return Property.js(name, TYPE.DATE_TIME, options); }

  static fromOptions(options: IProperty): Property { return Property.js(options.name!, options.type! as TYPE, options); }

  /** Registers the attached (dynamic) property for the specified type.
   * It is editable via the context panel, and gets saved into the view layout as well.
   * Property getter/setter typically uses Widget's "temp" property for storing the value. */
  static registerAttachedProperty(typeName: string, property: Property) {
    api.grok_Property_RegisterAttachedProperty(typeName, property.dart);
  }

  static propertyOptions:{[name in keyof IProperty]: IProperty & IPropertyMeta } = {
    'name': { name: 'name', type: TYPE.STRING, nullable: false },
    'type': { name: 'type', type: TYPE.STRING, nullable: false, description: 'Property data type, such as "int" or "string".' },
    'inputType': { name: 'inputType', type: TYPE.STRING, friendlyName: 'Input type', description: 'Property input type' },
    'nullable': { name: 'nullable', type: TYPE.BOOL, description: 'Whether an empty value is allowed. This is used by validators.' },
    'description': { name: 'description', type: TYPE.STRING, editor: InputType.TextArea, description: 'Property description' },
    'semType': { name: 'semType', type: TYPE.STRING, friendlyName: 'Semantic type', description: 'Semantic type' },
    'units': { name: 'units', type: TYPE.STRING, description: 'Units of measurement. See also: [postfix]' },
    'min': { name: 'min', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Minimum value. Applicable to numerical properties only' },
    'max': { name: 'max', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Maximum value. Applicable to numerical properties only' },
    'step': { name: 'step', applicableTo: TYPE.NUMERICAL, type: TYPE.FLOAT, description: 'Step to be used in a slider. Only applies to numerical properties.' },
    'showSlider': { name: 'showSlider', applicableTo: TYPE.NUMERICAL, type: TYPE.BOOL, description: 'Whether a slider appears next to the number input. Applies to numerical columns only.' },
    'showPlusMinus': { name: 'showPlusMinus', applicableTo: TYPE.NUMERICAL, type: TYPE.BOOL, description: 'Whether a plus/minus clicker appears next to the number input. Applies to numerical columns only.' },
    'choices': { name: 'choices', applicableTo: TYPE.STRING, type: TYPE.STRING_LIST, description: 'List of choices. Applicable to string properties only' },
    'initialValue': { name: 'initialValue', type: TYPE.STRING, description: 'Initial value used when initializing UI. See also {@link defaultValue}' },
    'defaultValue': { name: 'defaultValue', type: TYPE.OBJECT, description: 'Default value used for deserialization and cloning. See also {@link initialValue}.' },
    'editor': { name: 'editor', type: TYPE.STRING, description: 'Custom editor (such as slider or text area)' },
    'category': { name: 'category', type: TYPE.STRING, description: 'Corresponding category on the context panel' },
    'format': { name: 'format', type: TYPE.STRING, description: 'Value format, such as "0.00"' },
    'userEditable': { name: 'userEditable', type: TYPE.BOOL, description: 'Whether the property should be editable via the UI' },
    'validators': { name: 'validators', type: TYPE.STRING_LIST, description: 'List of validators. It can include [NAMED_VALIDATORS] as well as any pre-defined function names. Signature: validator(x: DG.Type): string | null. [null] indicates that the value is valid, [string] describes a validation error.' },
    'valueValidators': { name: 'valueValidators', type: TYPE.OBJECT, description: 'List of value validators (functions that take a value and return error message or null)' },
    'friendlyName': { name: 'friendlyName', type: TYPE.STRING, description: 'Custom field friendly name shown in [PropertyGrid]' },
    'fieldName': { name: 'fieldName', type: TYPE.STRING, description: 'Name of the corresponding JavaScript field. No need to specify it if it is the same as name.' },
    'tags': { name: 'tags', type: TYPE.MAP, description: 'Additional tags' },
    'options': { name: 'options', type: TYPE.MAP, description: 'Additional options.' },
  }
}


/** A dynamic property associated with the entity. */
export class EntityProperty extends Property {
  constructor(dart: any) {
    super(dart);
  }

  static create(name: string, type: string): EntityProperty {
    return toJs(api.grok_EntityProperty_Create(toDart(name), toDart(type)));
  }
}
