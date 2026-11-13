/**
 * Reaction drawing: RDKit's in both modes by the product owner's decision (crux has no `drawReaction` yet), the one
 * RDKit-only part of Chem's reactions. Kept apart from reaction-engine.ts, which runs reactions by engine; the grid's
 * reaction renderer is rendering/rdkit-reaction-renderer.ts.
 */
import {getRdKitModule} from '../chem-common-rdkit';
import {renderReactionToCanvas} from '../../rendering/rdkit-reaction-renderer';

/** Draws a reaction (or a route) on `canvas` for the reaction dialogs, the reaction browser and the editor. */
export function drawReactionByRdKit(canvas: HTMLCanvasElement, smarts: string, w?: number, h?: number): boolean {
  return renderReactionToCanvas(getRdKitModule(), canvas, smarts, w, h);
}
