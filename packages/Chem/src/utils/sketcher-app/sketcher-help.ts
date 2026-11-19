/**
 * The Sketcher app's help page, help/sketcher.md in Chem's files (the help pane reads it from the package, as it reads
 * help/reactions.md for the Reaction Sketcher): the app's view and its Export dialog each point at their section.
 */
import {chemWebRoot} from '../reactions/reactions-help';

/** The page, in the package's files. */
export const SKETCHER_HELP_FILE = 'help/sketcher.md';

/** The page's sections, as its headings' anchors. */
export type SketcherHelpSection = 'the-link' | 'export';

/** The page's address, at `section`. */
export function sketcherHelpUrl(section?: SketcherHelpSection): string {
  return `${chemWebRoot()}${SKETCHER_HELP_FILE}${section ? `#${section}` : ''}`;
}
