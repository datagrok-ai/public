/**
 * A reaction written into the Crux sketcher arranged (Arrange runs on demand and when a whole scheme is written in,
 * never by itself while drawing; without it the app's default route opened at 20% as one long row, hard to read). A
 * reaction or route text that carries no coordinates of its own is laid out in rows that break at its arrows
 * (crux-sketch's Arrange, API-070) when the Reaction Sketcher or the reaction cell editor writes it into the sketcher,
 * so it opens readable. A text with coordinates (an RXN or RD file, a CDXML page) keeps its own layout, as crux-sketch
 * keeps a pasted page's.
 */
import type {SetValueOptions, SketcherHandle} from '../../../vendor/crux-sketch/crux-sketch';
import {noReactionArrow} from './reaction-document';

/**
 * Whether the vendored sketcher's `setValue` takes `{arrange: true}`: the write itself arranged, with no undo step of
 * its own. Without it, `arrange()` follows the write, one edit the first Ctrl+Z takes back.
 */
export const SET_VALUE_ARRANGES = true;

/** Whether `text` is written arranged: it has a reaction arrow and no coordinates of its own. */
export function arrangesOnWrite(text: string): boolean {
  const t = text.trimStart();
  return t !== '' && !noReactionArrow(t) && !t.startsWith('$RXN') && !t.startsWith('$RDFILE') && !t.startsWith('<') &&
    !/\n\s*M {2}END/.test(t);
}

/** `text` into `sketch` (its `setValue`), arranged when {@link arrangesOnWrite} says so. Returns the sketcher's value
 * once arranged by an Arrange edit (as {@link arrangeOpened}), else null. */
export function writeArranged(sketch: SketcherHandle, text: string, options: SetValueOptions = {}): string | null {
  if (!arrangesOnWrite(text)) {
    sketch.setValue(text, options);
    return null;
  }
  if (SET_VALUE_ARRANGES) {
    sketch.setValue(text, {...options, arrange: true} as SetValueOptions);
    return null;
  }
  sketch.setValue(text, options);
  return sketch.lastError === null && arrangeNow(sketch) ? sketch.value : null;
}

/**
 * `sketch` opened on `text` (its first value, read): arranged when {@link arrangesOnWrite} says so, as no person's
 * edit, read as `format` when one is given (the cell editor's SMARTS column). Returns the sketcher's value once
 * arranged by an Arrange edit (no setValue option: Ctrl+Z takes it back, and a person's edits undone come back to it,
 * not to `text`), else null.
 */
export function arrangeOpened(sketch: SketcherHandle, text: string, format: 'smarts' | null = null): string | null {
  if (sketch.lastError !== null || !arrangesOnWrite(text))
    return null;
  if (SET_VALUE_ARRANGES) {
    sketch.setValue(text, {view: 'fit', arrange: true, ...(format ? {format} : {})} as SetValueOptions);
    return null;
  }
  return arrangeNow(sketch) ? sketch.value : null;
}

/** The sketcher's Arrange where the vendored sketcher has it; else nothing. */
export function arrangeNow(sketch: SketcherHandle): boolean {
  const arrange = (sketch as unknown as {arrange?: (width?: number) => boolean}).arrange;
  return typeof arrange === 'function' ? arrange.call(sketch) : false;
}
