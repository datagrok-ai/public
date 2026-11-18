/**
 * What a reaction written as text does not keep of the sketcher's drawing, so that it is named with a warning (a
 * reagent name drawn over the arrow, an atom label "HATU", came back as "*" in the app and in table cells, without a
 * word). Read from Crux Sketch's `warnings`, the rows of kind `format` for the notation the text is written in:
 * - `alias`: an atom drawn as a label (a name such as HATU). Reaction SMILES and SMARTS write it as a plain atom, "*"
 *   for a name, so calculations, search and mapping see no structure there; an RXN block keeps it.
 * - the drawing's own objects, which no reaction text writes (the objects of crux-sketch's reaction document): the
 *   texts, compound numbers, the arrows' styles, a failed reaction, resonance and curved arrows, brackets, cycles. A
 *   code Chem has no words for yet (a later object) is said in Crux Sketch's own words.
 * Molecule losses (stereo groups, queries, S-groups and the rest of crux-sketch's loss table) are left to the checks
 * that have rules for them. The words a person reads are in {@link DRAWING_LOSS_WORDS} and rules.ts.
 */
import type {SketchWarning} from '../../../vendor/crux-sketch/crux-sketch';
import {Flag, flag} from '../footguns/rules';
import {isCdxmlText} from '../../rendering/cdxml-reaction-cells';

/** The words for each drawing object a reaction text drops (Crux Sketch's loss codes), in the order they are said. */
export const DRAWING_LOSS_WORDS: {[code: string]: string} = {
  'texts': 'the texts (conditions and captions)',
  'compound-numbers': 'the compound numbers',
  'arrow-style': 'the arrows\' styles (each is written as a plain arrow)',
  'failed-reaction': 'the failed reaction (written as no reaction)',
  'resonance-arrows': 'the resonance arrows',
  'curved-arrows': 'the curved arrows',
  'brackets': 'the brackets (their repeat counts and charges)',
  'cycles': 'the catalytic cycles as drawn (each is written as its steps)',
};

/** Crux Sketch's molecule losses (its `SketchLossCode`s): not drawing objects, so not said here, but `alias`. */
const MOLECULE_LOSS_CODES = new Set(['r-label', 'dummy', 'dative', 'undefined-stereo', 'charges-rearranged',
  'unusual-valence', 'ambiguous-stereo', 'library', 'alias', 'r-label-numbers', 'attachment-points', 'stereo-groups',
  'chiral-flag', 'query', 'sgroup', 'hydrogen-bond', 'zero-order', 'atom-maps', 'stereo', 'v3000']);

/** The notation a reaction text is written in, as Crux Sketch's warnings name it: an RXN block, reaction SMARTS (SMARTS
 * syntax outside a CX part, as detectors.js reads it) or reaction SMILES. */
export type ReactionTextFormat = 'rxn' | 'smarts' | 'smiles' | 'cdxml';

/** The notation of reaction `text`. */
export function reactionTextFormat(text: string): ReactionTextFormat {
  if (text.trimStart().startsWith('$RXN'))
    return 'rxn';
  if (isCdxmlText(text))
    return 'cdxml';
  return /\[.?#\d|\$|&|;|,|!.?]|~/.test(text.replace(/\s+\|[^|]*\|\s*$/, '')) ? 'smarts' : 'smiles';
}

/** The words a flag names a notation with. */
const FORMAT_WORDS: {[format in ReactionTextFormat]: string} = {smiles: 'reaction SMILES', smarts: 'reaction SMARTS',
  rxn: 'an RXN block', cdxml: 'CDXML'};

/** What `text`, written from a drawing with `warnings`, does not keep of it. */
export interface DrawingLosses {
  /** The notation `text` is written in. */
  format: ReactionTextFormat;
  /** The atoms drawn as labels that the text writes as plain atoms (none in an RXN block, which keeps them). */
  labels: number;
  /** The drawing's objects the text drops, in words. */
  objects: string[];
}

/** What reaction `text`, written from a drawing whose Crux Sketch `warnings` these are, does not keep of it. The
 * drawing's objects are listed under SMILES and SMARTS alike (any reaction text drops them); a label's loss is read
 * under the text's own notation. */
export function drawingLosses(warnings: readonly SketchWarning[], text: string): DrawingLosses {
  const format = reactionTextFormat(text);
  if (format === 'cdxml') { // CDXML keeps labels and the drawing: what it does not keep is Crux Sketch's to say
    const said = warnings.filter((w) => w.kind === 'format' && w.format === 'cdxml').map((w) => w.message);
    return {format, labels: 0, objects: [...new Set(said)]};
  }
  const own = format === 'rxn' ? 'smiles' : format;
  const rows = warnings.filter((w) => w.kind === 'format' && w.format === own && typeof w.code === 'string');
  const aliasRows = rows.filter((w) => w.code === 'alias');
  // a row with no atoms listed still says one label at least
  const labels = Math.max(new Set(aliasRows.flatMap((w) => w.atoms)).size, aliasRows.length ? 1 : 0);
  const known = Object.keys(DRAWING_LOSS_WORDS);
  const codes = [...new Set(rows.map((w) => w.code as string).filter((c) => !MOLECULE_LOSS_CODES.has(c)))];
  // Chem's words in their order first, then any later object in Crux Sketch's words
  const objects = [...known.filter((c) => codes.includes(c)).map((c) => DRAWING_LOSS_WORDS[c]),
    ...codes.filter((c) => !known.includes(c)).map((c) => rows.find((w) => w.code === c)!.message)];
  return {format, labels: format === 'rxn' ? 0 : labels, objects};
}

/** The flags of `losses` (warns: OK still writes): a label written as a plain atom, then the drawing's
 * objects left out. */
export function drawingLossFlags(losses: DrawingLosses): Flag[] {
  const words = FORMAT_WORDS[losses.format];
  const flags: Flag[] = [];
  if (losses.labels > 0)
    flags.push(flag('labelAsStar', `${losses.labels === 1 ? 'An atom' : `${losses.labels} atoms`}|${words}`));
  if (losses.objects.length)
    flags.push(flag('drawingNotKept', `${words}: ${losses.objects.join(', ')}`));
  return flags;
}
