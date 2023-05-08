/// this file was generated automatically from d4 classes declarations
import { toDart } from "../wrappers";
let api = <any>window;

export class Inputs {
  static Int = 'Int';

  static BigInt = 'BigInt';

  static Float = 'Float';

  static QNum = 'QNum';

  static Slider = 'Slider';

  static Bool = 'Bool';

  static TextArea = 'TextArea';

  static Text = 'Text';

  static Date = 'Date';

  static Map = 'Map';

  static List = 'List';

  static Color = 'Color';

  static Column = 'Column';

  static Radio = 'Radio';

  static Molecule = 'Molecule';

  static UserGroupSelector = 'UserGroupSelector';

}
export function renderMultipleHistograms(g: CanvasRenderingContext2D, bounds: any, histograms: Array<Int32List>, options?: {categoryColumn?: any, colors?: Array<number>, tension?: number, normalize?: boolean, markerSize?: number, fill?: boolean}): any
  { return api.grok_renderMultipleHistograms(toDart(g), toDart(bounds), toDart(histograms), toDart(options?.categoryColumn), toDart(options?.colors), toDart(options?.tension), toDart(options?.normalize), toDart(options?.markerSize), toDart(options?.fill)); }

