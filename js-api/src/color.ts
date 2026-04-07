import type {Cell, Column} from "./dataframe";
import {MapProxy} from "./proxies";
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Color-related routines. */
export class Color {

    /** Returns a color associated with the specified cell as an ARGB-formatted integer.
     * To convert to html color, use {@link getCellColorHtml} or {@link toHtml}. */
    static getCellColor(cell: Cell): number {
        return api.grok_Color_FromCell(cell.dart);
    }

    /** Returns a string representation of the color associated with the specified cell.
     * For batch color manipulations, use {@link getCellColor}. */
    static getCellColorHtml(cell: Cell): string {
        return Color.toHtml(Color.getCellColor(cell));
    }

    /** Returns a color associated with the specified category within a column.
     * Returns ARGB-formatted integer. To convert to html color, use {@link toHtml}. */
    static getCategoryColor(column: Column, category: any): number {
        return api.grok_Color_FromCategory(column.dart, category);
    }

    /** Returns the Alpha component of the color represented as ARGB-formatted integer. */
    static a(c: number): number { return (c >> 24) & 0xFF; }

    /** Returns the Red component of the color represented as ARGB-formatted integer. */
    static r(c: number): number { return (c >> 16) & 0xFF; }

    /** Returns the Green component of the color represented as ARGB-formatted integer. */
    static g(c: number): number { return (c >> 8) & 0xFF; }

    /** Returns the Blue component of the color represented as ARGB-formatted integer. */
    static b(c: number): number { return c & 0xFF; }

    static argb(a: number, r: number, g: number, b: number) {
        return ((a << 24) | (r << 16) | (g << 8) | b) >>> 0;
    }

    /** Returns the color with the specified alpha component (0-255). */
    static setAlpha(color: number, alpha: number) {
        return Color.argb(alpha, Color.r(color), Color.g(color), Color.b(color));
    }

    /** Returns i-th categorical color (looping over the palette if needed) */
    static getCategoricalColor(i: number): number {
        return Color.categoricalPalette[i % Color.categoricalPalette.length];
    }

    /** Returns either black or white color, depending on which one would be most contrast to the specified [color]. */
    static getContrastColor(color: number): number {
        return api.grok_Color_GetContrastColor(color);
    }

    /** Converts ARGB-formatted integer color to a HTML-formatted string (such as `#ffffff`). See also {@link toRgb}. */
    static toHtml(color: number): string { return api.grok_Color_ToHtml(color); }

    /** Convert HTML-formatted string (such as `#ffffff`) to ARGB-fromatted integer color */
    static fromHtml(htmlColor: string): number { return api.grok_Color_FromHtml(htmlColor); }

    /** Converts ARGB-formatted integer color to a HTML-formatted string (such as `rbg(20, 46, 124)`). See also {@link toHtml. }*/
    static toRgb(color: number): string {
        return color === null ? '' : `rgb(${Color.r(color)},${Color.g(color)},${Color.b(color)})`;
    }

    /** For RDKit molecule substruct highlight */
    static hexToPercentRgb(hex: string): number[] | null {
        const result = hex.length === 7 ? /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex) :
            /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
        return result ? [
            parseInt(result[1], 16) / 256,
            parseInt(result[2], 16) / 256,
            parseInt(result[3], 16) / 256,
            result.length > 4 ? parseInt(result[4], 16) / 256 : 0.3
        ] : null;
    }


    /** Returns the standard palette of the categorical colors used across all visualizations in Datagrok. */
    static get categoricalPalette(): number[] {
        return api.grok_Color_CategoricalPalette();
    }

    /** Returns the map of existing palettes used in Datagrok. */
    static get categoricalPalettes(): {[key: string]: any} {
        return new MapProxy(api.grok_Color_GetCategoricalPalettes());
    }

    /** Returns the list of categorical color schemes used in Datagrok. */
    static get categoricalSchemes(): number[] {
        return api.grok_Color_CategoricalSchemes();
    }

    /** Returns the list of continuous color schemes for linear coloring used in Datagrok. */
    static get continuousSchemes(): number[] {
        return api.grok_Color_ContinuousSchemes();
    }


    static scaleColor(x: number, min: number, max: number, alpha?: number, colorScheme?: number[], belowMinColor?: number, aboveMaxColor?: number): number {
        return api.grok_Color_ScaleColor(x, min, max, alpha ? alpha : null, colorScheme ? colorScheme : null, belowMinColor ? belowMinColor : null, aboveMaxColor ? aboveMaxColor : null);
    }

    static highlight(color: number): number {
        return api.grok_Color_Highlight(color);
    }

    static darken(color: number, diff: number): number {
        return api.grok_Color_Darken(color, diff);
    }

    static  getRowColor(column: Column, row: number): number {
        return api.grok_Color_GetRowColor(column.dart, row);
    }

    static scale(x: number, min: number, max: number): number {
        return min === max ? min : (x - min) / (max - min);
    }

    static get gray(): number {
        return 0xFF808080;
    }

    static get lightLightGray(): number {
        return 0xFFF0F0F0;
    }

    static get lightGray(): number {
        return 0xFFD3D3D3;
    }

    static get darkGray(): number {
        return 0xFF838383;
    }

    static get blue(): number {
        return 0xFF0000FF;
    }

    static get green(): number {
        return 0xFF00FF00;
    }

    static get darkGreen(): number {
        return 0xFF006400;
    }

    static get black(): number {
        return 0xFF000000;
    }

    static get yellow(): number {
        return 0xFFFFFF00;
    }

    static get white(): number {
        return 0xFFFFFFFF;
    }

    static get red(): number {
        return 0xFFFF0000;
    }

    static get darkRed(): number {
        return 0xFF8b0000;
    }

    static get maroon(): number {
        return 0xFF800000;
    }

    static get olive(): number {
        return 0xFF808000;
    }

    static get orange(): number {
        return 0xFFFFA500;
    }

    static get darkOrange(): number {
        return 0xFFFF8C00;
    }

    static get lightBlue(): number {
        return 0xFFADD8E6;
    }

    static get darkBlue(): number {
        return 0xFF0000A0;
    }

    static get purple(): number {
        return 0xFF800080;
    }

    static get whitesmoke(): number {
        return 0xFFF5F5F5;
    }

    static get navy(): number {
        return 0xFF000080;
    }

    static get cyan(): number {
        return 0xFF00ffff;
    }

    static get filteredRows(): number {
        return 0xff1f77b4;
    }

    static get filteredOutRows(): number {
        return Color.lightLightGray;
    }

    static get selectedRows(): number {
        return Color.darkOrange;
    }

    static get missingValueRows(): number {
        return Color.filteredOutRows;
    }

    static get mouseOverRows(): number {
        return 0xFFAAAAAA;
    }

    static get currentRow(): number {
        return 0xFF38B738;
    }

    static get histogramBar(): number {
        return Color.filteredRows;
    }

    static get barChart(): number {
        return 0xFF24A221;
    }

    static get scatterPlotMarker(): number {
        return 0xFF40699c;
    }

    static get scatterPlotSelection(): number {
        return 0x80323232;
    }

    static get scatterPlotZoom(): number {
        return 0x80626200;
    }

    static get areaSelection(): number {
        return Color.lightBlue;
    }

    static get rowSelection(): number {
        return 0x60dcdca0;
    }

    static get colSelection(): number {
        return 0x60dcdca0;
    }

    static get areaZoom(): number {
        return 0x80323232;
    }

    static get gridWarningBackground(): number {
        return 0xFFFFB9A7;
    }

    static get success(): number {
        return 0xFF3cb173;
    }

    static get failure(): number {
        return 0xFFeb6767;
    }
}
