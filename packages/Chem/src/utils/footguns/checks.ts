/**
 * Foot-gun checks as Chem runs them (the rules, their levels and words: rules.ts): in Crux mode only (RDKit mode stays
 * as it was, so every check here answers nothing there), with Crux's module, for a value of a column (the context
 * panel's Checks section, widgets/checks-panel.ts) or a value an editor is about to write (the dialog footer's middle
 * slot, footer-slot.ts). The reaction cell editor's checks run in both modes, as it is the Crux editor whatever the
 * engine ({@link reactionEditorModule}).
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {ChemicalEngine, engineInEffect} from '../../engine/chem-engine';
import {cruxLoadFailed, cruxModuleSync, noteCruxError} from '../../engine/crux-module';
import {tooLargeToDraw} from '../../rendering/crux-reaction-drawing';
import {setEnumeratorErrorHook} from './enumerator';
import type {CruxReadFailure} from '../../crux/crux-sketcher';
import type {SketchWarning} from '../../../vendor/crux-sketch/crux-sketch';
import {drawingLossFlags, drawingLosses} from '../reactions/reaction-losses';
import {Flag, MoleculeTarget, ReactionTarget, flag, isMolBlockText, mappingFlags, moleculeFlags,
  reactantCountFlags, reactionEditFlags, reactionFlags, readsAsTemplate, setCruxGuards, sketcherWriteFlags,
  templateMappingFlags} from './rules';

// every error of a Crux call the checks make is heard (a trap the module could not recover from turns Crux off for the
// session), and no Crux call is made while it is off
setCruxGuards({noteError: noteCruxError, failed: cruxLoadFailed});
setEnumeratorErrorHook(noteCruxError);

export type {Flag, Level} from './rules';
export {flagsText, worstLevel, LEVEL_ORDER} from './rules';

/** Crux's module while the checks run (Crux mode, its module ready and not failed), else null: RDKit mode is never
 * flagged, and a session whose Crux trapped makes no more Crux calls. */
export function checksModule(): RDModule | null {
  return engineInEffect() === ChemicalEngine.Crux && !cruxLoadFailed() ? cruxModuleSync() : null;
}

/** Crux's module for the reaction cell editor's checks, whatever the Chemical Engine (the reaction editor is the Crux
 * editor in both modes and behaves the same, its footer's checks too; the editor loads
 * Crux in RDKit mode, crux-reaction-drawing.ts cruxForReactions): null until it is loaded, and once it has failed. The
 * molecule editor's checks, the Checks section and every other path keep {@link checksModule}'s Crux mode only. */
export function reactionEditorModule(): RDModule | null {
  return cruxLoadFailed() ? null : cruxModuleSync();
}

/** Whether a value is too large to check synchronously (the limit for drawing: a molecule over 2,500
 * characters, a whole over 10,000): the checks say nothing of it, rather than freeze the tab. */
function tooLarge(value: string): boolean {
  return tooLargeToDraw(value);
}

/** A molecule column's target: its units (a column without them holds SMILES). */
export function moleculeTarget(col: DG.Column): MoleculeTarget {
  return {units: col.meta.units || DG.chem.Notation.Smiles};
}

/** A reaction column's target: a template column holds reaction SMARTS (units `smarts`, as detectors.js writes them);
 * its units kept, for the note on a template in a column taken for something else. */
export function reactionTarget(col: DG.Column): ReactionTarget {
  return {template: col.meta.units === 'smarts', units: col.meta.units || undefined};
}

/** The flags of a molecule sketcher's OK into a column of `target` (rules.ts sketcherWriteFlags): `written` what OK
 * would write, `smilesQuery` non-null when the sketcher holds a query no SMILES writes, `drawing` its molblock,
 * `previous` the cell's value as the editor opened (what the edit added is flagged, what the cell held is not). */
export function checkSketcherWrite(written: string, target: MoleculeTarget, smilesQuery: string | null,
  drawing: string, previous = ''): Flag[] {
  const module = checksModule();
  return module && !tooLarge(written) ? sketcherWriteFlags(module, written, target, smilesQuery, drawing, previous) :
    [];
}

/** The flags of what a molecule sketcher's last read left out (crux-sketch's `readWarnings` codes): a reaction's arrow
 * (MOL-17). */
export function checkReadLosses(codes: readonly string[]): Flag[] {
  return checksModule() && codes.includes('arrow') ? [flag('arrowLeftOut')] : [];
}

/** The notations by the words a flag names them with. */
const NOTATION_WORDS: {[notation: string]: string} = {smiles: 'SMILES', molblock: 'a V2000 molblock',
  molblockV3000: 'a V3000 molblock', smarts: 'SMARTS'};

/** The block of a Crux getter that threw for what OK would write (the engine-poison foot-gun): the engine stopped, or
 * the drawing cannot be written in the column's notation (named; the sketcher's own message, which names the engine,
 * goes to the console, as the footer names no engine). */
export function checkReadFailure(failure: CruxReadFailure | null): Flag[] {
  if (!checksModule() || failure === null)
    return [];
  if (!failure.engine)
    console.warn(`Chem | the drawing could not be written as ${failure.notation}: ${failure.message}`);
  return [failure.engine ? flag('engineStopped') :
    flag('notWritable', NOTATION_WORDS[failure.notation] ?? failure.notation)];
}

/** The flags of reaction `text` in a column of `target` (a reaction editor's template: `{template: true}`). */
export function checkReaction(text: string, target: ReactionTarget): Flag[] {
  const module = checksModule();
  return module && !tooLarge(text) ? reactionFlags(module, text, target) : [];
}

/** The flags of the reaction cell editor's OK (rules.ts reactionEditFlags): `original` the cell's text as it opened,
 * `text` what OK would write, `edited` whether a person changed it. */
export function checkReactionEdit(original: string, text: string, target: ReactionTarget, edited: boolean): Flag[] {
  const module = reactionEditorModule(); // the reaction editor's: in both modes
  return module && !tooLarge(text) && !tooLarge(original) ?
    reactionEditFlags(module, original, text, target, edited) : [];
}

/** The flags of what reaction `text`, written from a drawing with Crux Sketch's `warnings`, does not keep of it (a
 * label written as "*", the drawing's texts and arrow styles; reactions/reaction-losses.ts). */
export function checkDrawingLosses(warnings: readonly SketchWarning[], text: string): Flag[] {
  // the reaction editor's: in both modes
  return reactionEditorModule() && text.trim() !== '' ? drawingLossFlags(drawingLosses(warnings, text)) : [];
}

/** The flags of running template `smarts` in a Transform dialog that gives it `reactants` reactants: its reactant
 * count, then its template checks (agents over the arrow left out: a Transform runs them, only enumeration skips). */
export function checkTransform(smarts: string, reactants: number): Flag[] {
  return checksModule() ? [...reactantCountFlags(smarts, reactants),
    ...checkReaction(smarts, {template: true}).filter((f) => f.rule !== 'agentsOverArrow')] : [];
}

/** The flags of an atom mapping at `confidence`, and of its reaction's products holding `missing`, heavy atoms no
 * reactant or reagent drawn holds ("C 2, O 1"; the mapper maps them onto something else, however sure it is),
 * when some are. */
export function checkMapping(confidence: number, missing: string | null = null): Flag[] {
  return checksModule() ? [...mappingFlags(confidence), ...(missing ? [flag('unbalancedMapping', missing)] : [])] : [];
}

/** Whether reaction value `text` is a template (MAP-02): its column's units are `smarts` (the current cell's, when
 * the context panel shows it), or it is one by what it holds (rules.ts readsAsTemplate). */
export function isTemplateValue(module: RDModule, text: string): boolean {
  const o = grok.shell.o;
  if (o instanceof DG.SemanticValue && o.value === text && o.cell?.column?.semType === DG.SEMTYPE.CHEMICAL_REACTION &&
    reactionTarget(o.cell.column).template)
    return true;
  return readsAsTemplate(module, text);
}

/** The flags of Map atoms on reaction value `text` when it is a template (MAP-02), else none. */
export function checkTemplateMapping(text: string): Flag[] {
  const module = checksModule();
  return module && !tooLarge(text) && isTemplateValue(module, text) ? templateMappingFlags() : [];
}

/** The flags of a cell's value (the context panel's Checks section): a molecule or a reaction, by its column. */
export function checkCell(cell: DG.Cell): Flag[] {
  const module = checksModule();
  const value = cell?.value;
  const col = cell?.column;
  if (!module || typeof value !== 'string' || !value || !col || tooLarge(value))
    return [];
  if (col.semType === DG.SEMTYPE.MOLECULE)
    return moleculeFlags(module, value, moleculeTarget(col));
  // a CDXML page (units `cdxml`) is no reaction text: Crux Sketch reads it, and its checks are its own
  if (col.semType === DG.SEMTYPE.CHEMICAL_REACTION)
    return col.meta.units === 'cdxml' ? [] : reactionFlags(module, value, reactionTarget(col));
  return [];
}

/** The flags of a semantic value the context panel shows: its cell's, by its column; a value with no cell (or one the
 * panel's condition sees only as a string, `value` of `semType`) as a SMILES or molblock, or a reaction that is no
 * template. */
export function checkValue(value: string, semType: string, cell?: DG.Cell | null): Flag[] {
  const module = checksModule();
  if (!module || !value || tooLarge(value))
    return [];
  if (!cell) {
    // the context panel shows the current object: its cell, when it is the one whose value this is
    const o = grok.shell.o;
    if (o instanceof DG.SemanticValue && o.value === value && o.cell?.column?.semType === semType)
      cell = o.cell;
  }
  if (cell?.column)
    return checkCell(cell);
  if (semType === DG.SEMTYPE.MOLECULE) {
    return moleculeFlags(module, value,
      {units: isMolBlockText(value) ? DG.chem.Notation.MolBlock : DG.chem.Notation.Smiles});
  }
  if (semType === DG.SEMTYPE.CHEMICAL_REACTION)
    return /<CDXML[\s>]/.test(value.slice(0, 4096)) ? [] : reactionFlags(module, value, {template: false});
  return [];
}
