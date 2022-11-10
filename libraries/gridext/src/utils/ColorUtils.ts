export class ColorUtils {

  static a(c: number): number {
    return (c >> 24) & 0xFF;
  }

  /** Returns the Red component of the color represented as ARGB-formatted integer. */
  static r(c: number): number {
    return (c >> 16) & 0xFF;
  }

  /** Returns the Green component of the color represented as ARGB-formatted integer. */
  static g(c: number): number {
    return (c >> 8) & 0xFF;
  }

  /** Returns the Blue component of the color represented as ARGB-formatted integer. */
  static b(c: number): number {
    return c & 0xFF;
  }

  static toRgb(color: number): string {
    return color === null ? '' : `rgb(${ColorUtils.r(color)},${ColorUtils.g(color)},${ColorUtils.b(color)})`;
  }

  static get currentRow(): number {
    return 0xFF38B738;
  }

  static get colSelection(): number {
    return 0x60dcdca0;
  }

  static get sortArrow(): number {
    return 0xFF0000A0;
  }
}
