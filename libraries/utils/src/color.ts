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

export function toRgba(color: DG.ColorType): string {
  return `rgba(${DG.Color.r(color)}, ${DG.Color.g(color)}, ${DG.Color.b(color)}, ${DG.Color.a(color) / 255})`;
}
