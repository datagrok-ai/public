import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** Set alpha over {@link color} */
export function setAlpha(color: DG.ColorType, alpha: number): number {
  return 0 |
    Math.round(DG.Color.a(color) * alpha) << 24 |
    DG.Color.r(color) << 16 |
    DG.Color.g(color) << 8 |
    DG.Color.b(color);
}

export function intToRgb(value: DG.ColorType): string {
  return `rgb(${DG.Color.r(value)}, ${DG.Color.g(value)}, ${DG.Color.b(value)})`;
}

export function intToRgba(value: DG.ColorType): string {
  return `rgba(${DG.Color.r(value)}, ${DG.Color.g(value)}, ${DG.Color.b(value)}, ${DG.Color.a(value) / 255})`;
}

export function intToHtml(value: DG.ColorType): string {
  return `#${(value & 0xFFFFFF).toString(16).padStart(6, '0')}`;
}

export function intToHtmlA(value: DG.ColorType): string {
  const alpha = DG.Color.a(value);
  return `#${(value & 0xFFFFFF).toString(16).padStart(6, '0')}` +
    alpha.toString(16).padStart(2, '0');
}

/** Converts color Hex string ('#000bff') to int number */
export function htmlToInt(value: string): number {
  return DG.Color.fromHtml(value) & 0xFFFFFF;
}

/** Converts color Hex string ('#000bff22' , alpha=0x22) to int number with alpha */
export function htmlToIntA(value: string): number {
  const m: RegExpExecArray = /#([0-9A-Fa-f]{6})([0-9A-Fa-f]{2})/.exec(value)!;
  return parseInt(m[1], 16) | (parseInt(m[2], 16) << 24);
}
