import type { Density } from './types.js';
export type { Density };
/** The density the bars are drawn at. */
export type BarDensity = Exclude<Density, 'auto'>;
/** Under this width or this height, in CSS px of the sketcher's box, an `'auto'` sketcher's bars are compact (B1). */
export declare const COMPACT_BELOW: Readonly<{
    width: 720;
    height: 600;
}>;
/** A box, in CSS px. */
export interface BoxSize {
    readonly width: number;
    readonly height: number;
}
/**
 * The option as written, by HTML's rule for enumerated attributes: `'regular'` or `'compact'` only when set to that
 * keyword (in any case, spaces around it ignored), and `'auto'` for anything else, nothing set included.
 */
export declare function toDensity(value: unknown): Density;
/**
 * The density the bars are drawn at: the option's, unless it is `'auto'`; then compact in a box under 720 px wide or
 * under 600 px tall, else regular. A box not measured yet, or with no size (a hidden sketcher), is drawn regular.
 */
export declare function densityOf(option: Density, box: BoxSize | null): BarDensity;
