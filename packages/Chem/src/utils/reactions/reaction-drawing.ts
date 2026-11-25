/**
 * Reaction previews: the reaction dialogs' (Transformation, Two-Component and their views, full screen included), the
 * reaction browser's cards and the reaction editor's. They are drawn by the engine in effect, as the reaction cells
 * are (rendering/rdkit-reaction-renderer.ts): in Crux mode by Crux (rendering/crux-reaction-drawing.ts), with RDKit
 * behind it for a template Crux reads no reaction from (the fallback `reaction preview`, engine/fallback.ts); in RDKit
 * mode by RDKit, as before. Kept apart from reaction-engine.ts, which runs reactions by engine.
 */
import {getRdKitModule} from '../chem-common-rdkit';
import {renderReactionToCanvas} from '../../rendering/rdkit-reaction-renderer';
import {chemDrawOptions, cruxDrawsReactions, cruxFittedPicture, cruxReadsReaction, drawReactionByCrux,
  tooLargeToDraw} from '../../rendering/crux-reaction-drawing';
import {FULL_SCREEN_CENTRE, FULL_SCREEN_FILL, FULL_SCREEN_LINE_SHARE, FULL_SCREEN_MAX_BOND}
  from '../../widgets/reaction-panel';
import {cruxOrRdKit} from '../../engine/fallback';

/** The fallback's call name: a preview RDKit drew for Crux, logged once a session as
 * "Crux fallback: reaction preview". */
export const PREVIEW_CALL = 'reaction preview';

/** How a preview is drawn besides its size. */
export interface PreviewOptions {
  /** A full-screen view: in Crux mode the reaction is fitted to the box, as the Reaction panel's full-screen view
   * fits it (its bond grows up to {@link FULL_SCREEN_MAX_BOND} px, its lines a share of the bond); else Crux's own
   * bond, a reaction too large for the box shrunk to fit it, as in the cells. RDKit fills the box either way. */
  fullScreen?: boolean;
}

/** Tests: the texts Crux is taken to read no reaction from (see {@link setCruxUnreadablePreviewsForTests}). */
let unreadableForTests: Set<string> | null = null;

/** Tests: Crux reads no reaction from these texts in the previews (null: none); RDKit draws them, as the fallback. */
export function setCruxUnreadablePreviewsForTests(texts: string[] | null): void {
  unreadableForTests = texts ? new Set(texts) : null;
}

/**
 * Draws `reaction` (a reaction or a route) on the whole of `canvas`, `w` x `h` device px (the canvas's own size by
 * default) at the page's device pixel ratio, cleared first: by Crux in Crux mode, by RDKit in RDKit mode, as the
 * reaction cells. False when RDKit reads no reaction from it (RDKit mode, or Crux's fallback).
 */
export function drawReactionPreview(canvas: HTMLCanvasElement, reaction: string, w: number = canvas.width,
  h: number = canvas.height, options: PreviewOptions = {}): boolean {
  const rdkit = () => renderReactionToCanvas(getRdKitModule(), canvas, reaction, w, h);
  if (!cruxDrawsReactions())
    return rdkit();
  return cruxOrRdKit(PREVIEW_CALL, () => {
    // a reaction too large to draw is said to be, as in the cells (drawReactionByCrux)
    if (!tooLargeToDraw(reaction) && (!cruxReadsReaction(reaction) || unreadableForTests?.has(reaction)))
      throw new Error(`Crux reads no reaction from ${reaction}`);
    const g = canvas.getContext('2d')!;
    g.save();
    g.setTransform(1, 0, 0, 1, 0, 0);
    g.clearRect(0, 0, w, h);
    g.restore();
    const r = window.devicePixelRatio;
    const [width, height] = [w / r, h / r];
    const drawOptions = chemDrawOptions();
    if (!options.fullScreen) {
      drawReactionByCrux(g, reaction, 0, 0, width, height, {drawOptions});
      return true;
    }
    const fit = cruxFittedPicture(reaction, width, height, {drawOptions, centre: FULL_SCREEN_CENTRE},
      FULL_SCREEN_MAX_BOND, FULL_SCREEN_FILL);
    const lineWidth = fit ? Math.max(Number(drawOptions.bondLineWidth ?? 1),
      Math.round(fit.bondLength * FULL_SCREEN_LINE_SHARE * 10) / 10) : undefined;
    drawReactionByCrux(g, reaction, 0, 0, width, height, {
      drawOptions: lineWidth === undefined ? drawOptions : {...drawOptions, bondLineWidth: lineWidth},
      bondLength: fit?.bondLength, columnWidth: fit?.columnWidth, centre: FULL_SCREEN_CENTRE});
    return true;
  }, rdkit);
}
