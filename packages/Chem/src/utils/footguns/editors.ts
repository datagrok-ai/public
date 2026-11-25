/**
 * The foot-gun flags of Chem's editor dialogs, in their footer's middle slot (footer-slot.ts): the molecule cell
 * editor (package.ts editMoleculeCell), the reaction cell editor (reactions/reaction-cell-editor.ts), New/Edit Reaction
 * (reactions/reaction-editor.ts) and the Transform dialogs (reactions/ui.ts). Crux mode only, but for the reaction cell
 * editor, which is Crux's in both modes: in RDKit mode the other dialogs get no slot and nothing changes. Also the
 * molecule cell editor's "was it edited" rule (MOL-07).
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {CruxSketcher} from '../../crux/crux-sketcher';
import {sameMolecule} from '../chem-common-rdkit';
import {noteCruxError} from '../../engine/crux-module';
import {FooterSlot, footerSlot} from './footer-slot';
import {checkDrawingLosses, checkReactionEdit, checkReadFailure, checkReadLosses, checkSketcherWrite, checksModule,
  moleculeTarget, reactionTarget} from './checks';
import type {ReactionCellEditor} from '../reactions/reaction-cell-editor';
import type {CruxReadFailure} from '../../crux/crux-sketcher';


/**
 * The molecule cell editor's slot: after each change of the sketcher (and once it is ready) the flags of what OK would
 * write into `col` (rules.ts sketcherWriteFlags): a block disables OK. `previous` is the cell's value as the editor
 * opened. Null in RDKit mode (no slot).
 */
export function watchMoleculeEditor(dialog: DG.Dialog, sketcher: DG.chem.Sketcher, col: DG.Column, previous: string):
  FooterSlot | null {
  if (!checksModule())
    return null;
  const slot = footerSlot(dialog);
  const target = moleculeTarget(col);
  const smiles = target.units === DG.chem.Notation.Smiles;
  const update = () => {
    // a query Crux's SMILES can't hold is asked first: reading that SMILES would raise the platform's warning
    const impl = sketcher.sketcher;
    const query = smiles && impl instanceof CruxSketcher && impl.smilesIsQuery;
    const written = query ? '' : smiles ? sketcher.getSmiles() : sketcher.getMolFile();
    const drawing = sketcher.getMolFile();
    // a getter that threw for what OK would write (the engine stopped, or the drawing can't be written so: it answered
    // ''): a block, and OK writes nothing (the engine-poison foot-gun)
    const failed = checkReadFailure(moleculeEditorReadFailure(sketcher, target.units));
    if (failed.length) {
      slot.show(failed);
      return;
    }
    // what a paste left out (a ChemDraw reaction's arrow: MOL-17), then what OK would write
    const losses = impl instanceof CruxSketcher ? checkReadLosses(impl.readLosses) : [];
    slot.show([...losses, ...checkSketcherWrite(written, target, query ? '' : null, drawing, previous)]);
  };
  const changes = sketcher.onChanged.subscribe(update);
  sketcher.sketcherReady().then(update);
  const closed = dialog.onClose.subscribe(() => {
    changes.unsubscribe();
    closed.unsubscribe();
  });
  update();
  return slot;
}

/**
 * Why the molecule cell editor's OK must write nothing (the engine-poison foot-gun): Crux's getter for the column's
 * notation threw on its last read (the engine stopped, or the drawing cannot be written so), so the value the host
 * gives is '' and no edit. Null otherwise, in RDKit mode (which stays as it was) and for another sketcher.
 */
export function moleculeEditorReadFailure(sketcher: DG.chem.Sketcher, units: string | null | undefined):
  CruxReadFailure | null {
  const impl = sketcher.sketcher;
  if (!checksModule() || !(impl instanceof CruxSketcher))
    return null;
  return impl.readFailure(units === DG.chem.Notation.Smiles ? ['smiles'] : ['molblock', 'molblockV3000']);
}

/**
 * Whether the molecule cell editor's new value is the cell's own, so that OK writes nothing (a geometry-only edit keeps
 * the cell as it was). Crux mode compares canonical CXSMILES without coordinates, so an edit of only what a canonical
 * SMILES drops (enhanced stereo, an attachment point, an S-group) is written (MOL-07); RDKit mode compares canonical
 * SMILES, as it always has.
 */
export function unchangedByEdit(before: string, after: string): boolean {
  const module = checksModule();
  if (!module)
    return sameMolecule(before, after);
  const cx = (m: string): string | null => {
    let mol: RDMol | null = null;
    try {
      mol = module.get_mol(m);
      return mol ? mol.get_cxsmiles(JSON.stringify({CX_ALL_BUT_COORDS: true})) : null;
    } catch (e) {
      noteCruxError(e); // a trap the module could not recover from turns Crux off
      return null;
    } finally {
      mol?.delete();
    }
  };
  const a = cx(before);
  // S-groups too: Crux keeps them in Chem's stereo setting since chem-wasm 732af7a03 (before it, it read none, and
  // a text comparison of the S-groups stood in)
  return a !== null && a === cx(after);
}

/** The reaction cell editor's slot: once the editor has read the cell, what it read wrongly (a value
 * it takes for a molecule: RXN-15); after each change, the flags of the text OK would write into the cell's column (no
 * arrow left, the template checks, query features lost, an RXN block written otherwise: rules.ts reactionEditFlags);
 * a reaction Crux Sketch did not open (`refused`: RXN-06) a block throughout. In both modes (the reaction editor is
 * the Crux editor whatever the Chemical Engine; checks.ts reactionEditorModule): its flags show
 * once Crux is loaded, which the editor waits for before it reads the cell. */
export function watchReactionEditor(editor: ReactionCellEditor): FooterSlot {
  const slot = footerSlot(editor.dialog);
  const target = reactionTarget(editor.cell.cell.column);
  // RXN-06: a reaction Crux Sketch did not open (a convergent route) leaves OK off whatever happens after, as OK
  // writes nothing then (the editor's `write`); the editor says why once, in place of its sketcher, so the slot does
  // not say it again (it would be said twice, above the drawing and in the footer): the cell's own
  // value is what stays, so its flags are the unedited ones. Otherwise an edit's text names what it does not keep of
  // the drawing (a label written as "*", the texts); a CDXML cell's page is no reaction text to check, only
  // what CDXML does not keep of the drawing
  const update = (text: string) => editor.refused !== null ?
    slot.show(checkReactionEdit(editor.original, editor.original, target, false), true) :
    slot.show([...(editor.holdsCdxml ? [] : checkReactionEdit(editor.original, text, target, editor.edited)),
      ...(editor.edited ? checkDrawingLosses(editor.warnings, text) : [])]);
  editor.onChange.subscribe(update);
  editor.ready.then(() => {
    if (editor.dialog.root.isConnected)
      update(editor.value);
  });
  return slot;
}
