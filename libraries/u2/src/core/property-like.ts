/** The property metadata a form is generated from, and the shape the spec registry describes a
 * component's props with. `DG.Property` IS a `PropertyLike`; so is any descriptor that names the
 * field and can read and write it (`type` is read when `propertyType` is absent, as on a
 * `DG.IProperty` literal). A property with no `set` renders read-only. */
export interface PropertyLike {
  name: string;
  caption?: string | null;
  friendlyName?: string | null;
  propertyType?: string | null;
  type?: string | null;
  semType?: string | null;
  description?: string | null;
  /** Groups the property in a property panel; ungrouped properties render together. */
  category?: string | null;
  choices?: string[] | null;
  /** Platform `InputType` ('Radio', 'Color', 'Tags'…): the editor named here wins over every
   * other hint, as it does in `input_base.dart:668`. */
  inputType?: string | null;
  /** Platform editor hint ('textarea', 'password', 'switch', 'slider'), consulted next, and only
   * for the types Dart honors it for (`input_base.dart:702-728`). */
  editor?: string | null;
  nullable?: boolean;
  min?: number | null;
  max?: number | null;
  step?: number | null;
  /** A platform format string ('#0.00'), applied through `DG.format` where the platform is loaded. */
  format?: string | null;
  units?: string | null;
  showSlider?: boolean | null;
  showPlusMinus?: boolean | null;
  get?: (source: any) => any;
  set?: (source: any, value: any) => void;
}
