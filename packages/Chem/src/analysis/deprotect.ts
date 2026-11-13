/* eslint-disable max-len */
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Observable, Subject} from 'rxjs';
import {hasNewLines} from '../utils/chem-common';
import {MAX_SMILES_LENGTH} from '../utils/chem-constants';
import {EngineReaction} from '../utils/reactions/reaction-engine';
import {byEngineModule} from '../engine/chem-engine';
import {deprotectOne} from '../crux/jobs/reactions';
import {getCruxService} from '../crux/crux-searches';

/** A Deprotect job's time per molecule on the Crux workers (ms; CruxService.jobInChunks). */
const DEPROTECT_MS_PER_MOLECULE = 20;
/** Tests (benchmarks): false runs the column on the main thread, as before. */
let deprotectOnWorkers = true;

/** Tests (benchmarks): Deprotect's column on the Crux workers (true, the default) or on the main thread. */
export function setDeprotectOnWorkersForTests(on: boolean): void {
  deprotectOnWorkers = on;
}

/** The molecules' results on the Crux workers (crux/jobs/reactions.ts), null where the main thread is to react one
 * (an empty cell, a molecule Crux threw on); null for all when the workers fail. */
async function deprotectByWorkers(molecules: (string | null)[], smirks: string): Promise<(string | null)[] | null> {
  try {
    const parts = await (await getCruxService()).jobInChunks<(string | null)[]>('deprotect', molecules,
      DEPROTECT_MS_PER_MOLECULE, [smirks]);
    return parts.flat();
  } catch (e) {
    console.warn(`Chem | Deprotect on the Crux workers failed, the column is deprotected here: ${e}`);
    return null;
  }
}

function correctRGroups(smiles: string): string {
  const elementRGroupRegex = /\[R[1-9]\]/g;
  // replace all [R1] with [*:1]
  let correctedSmiles = smiles.replaceAll(elementRGroupRegex, (match) => {
    const rGroupNum = match[2];
    return `[*:${rGroupNum}]`;
  });

  // in some scenarios, rgroups can be written as [2*]
  const elementRGroupRegex2 = /\[\d\*\]/g;
  correctedSmiles = correctedSmiles.replaceAll(elementRGroupRegex2, (match) => {
    const rGroupNum = match[1];
    return `[*:${rGroupNum}]`;
  });

  // in some other scenarios, rgroups can be written as [1*:1] or [1*:0]
  const elementRGroupRegex3 = /\[\d\*\:\d\]/g;
  return correctedSmiles.replaceAll(elementRGroupRegex3, (match) => {
    const rGroupNum = match[1];
    return `[*:${rGroupNum}]`;
  });
}

function radicalToRgroupForm(smiles: string): string {
  // here we are sure that all marked carbons are in the form [C:1], [C:2], ..., [C:9]
  // replace all [C:1] with [*:1]
  return smiles.replaceAll(/\[C\:\d\]/g, (match) => {
    const rGroupNum = match[3];
    return `[*:${rGroupNum}]`;
  });
}

export const defaultFragment = 'O=C([*:1])OCC1c2ccccc2-c2ccccc21';
export class DeprotectEditor extends DG.FuncCallEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  columnInput!: DG.InputBase<DG.Column | null>;
  fragmentInput: DG.InputBase<string>;
  inputChangedSubject: Subject<any> = new Subject<any>();
  constructor(private funcCall: DG.FuncCall) {
    const root = ui.divV([]);
    super(root);
    this.tableInput = ui.input.table('Table', {value: funcCall.inputs['table'] ?? grok.shell.t, nullable: false, onValueChanged: () => {
      this.updateColumnInput();
      funcCall.inputs['table'] = this.tableInput.value;
      this.inputChangedSubject.next();
    }, tooltipText: 'Input data frame containing molecule column'});
    this.fragmentInput = ui.input.molecule('Fragment', {value: radicalToRgroupForm(funcCall.inputs['fragment'] ?? defaultFragment), nullable: false, onValueChanged: () => {
      funcCall.inputs['fragment'] = correctRGroups(grok.chem.convert(this.fragmentInput.value ?? '', grok.chem.Notation.Unknown, grok.chem.Notation.Smiles));
      this.inputChangedSubject.next();
    }, tooltipText: 'Smiles of the protecting group fragment to be removed. R1 group should be used to mark the attachment point.'});
    root.appendChild(this.tableInput.root);
    root.appendChild(this.fragmentInput.root);
    this.updateColumnInput();
    this.fragmentInput.addValidator(() => this.validate());
    this.tableInput.addValidator(() => this.validate());
    // trigger changes to init the funccall
    this.tableInput.fireChanged();
    this.fragmentInput.fireChanged();
    this.columnInput.fireChanged();
  }

  updateColumnInput() {
    const table = this.tableInput.value;
    if (table == null)
      return;
    const molColumn = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.columnInput?.root?.remove();
    this.columnInput = ui.input.column('Molecules', {
      table: table,
      nullable: false, filter: (col) => col.semType === DG.SEMTYPE.MOLECULE, value: molColumn ?? undefined, onValueChanged: () => {
        this.funcCall.inputs['molecules'] = this.columnInput.value;
        this.inputChangedSubject.next();
      }, tooltipText: 'Column with molecules to deprotect'});
    this.root.insertBefore(this.columnInput.root, this.fragmentInput.root);
    this.columnInput.addValidator(() => this.validate());
  }

  validate = () => {
    if (!this.tableInput.value || !this.columnInput.value || !this.fragmentInput.value)
      return 'All inputs must be set';
    const correctedFragment = correctRGroups(grok.chem.convert(this.fragmentInput.value ?? '', grok.chem.Notation.Unknown, grok.chem.Notation.Smiles));
    if (!correctedFragment)
      return 'Fragment must be a valid molecule';
    if (correctedFragment.indexOf('[*:1]') === -1)
      return 'Fragment must contain R1-group';
    return null;
  };

  get isValid(): boolean {
    // validate inputs
    return this.tableInput.validate() && this.columnInput.validate() && this.fragmentInput.validate();
  }

  inputFor(propertyName: string): DG.InputBase {
    switch (propertyName) {
    case 'table':
      return this.tableInput;
    case 'molecules':
      return this.columnInput;
    case 'fragment':
      return this.fragmentInput;
    default:
      throw new Error(`Unknown property name: ${propertyName}`);
    }
  }

  get onInputChanged(): Observable<any> {
    return this.inputChangedSubject;
  }
}

/** In Crux mode the reaction runs on the Crux module, the column on the Crux workers, and on RDKit's (`rdModule`) per
 * molecule Crux cannot do. */
export async function addDeprotectedColumn(table: DG.DataFrame, molColumn: DG.Column, fragment: string, rdModule: RDModule): Promise<DG.Column | null> {
  try {
    // "Chemical engine": the fragment by the main thread's module (Crux's in Crux mode, which runs the reaction), RDKit's
    // where Crux reads none
    const hFragmentSmiles = byEngineModule('deprotection fragment', rdModule, (module) => {
      const fragmentMol = module.get_mol(fragment);
      if (!fragmentMol)
        return null;
      try {
        // we need to add explicit hydrogens to the fragment to make sure the reaction works correctly
        fragmentMol.add_hs_in_place();
        return fragmentMol.get_smiles();
      } finally {
        fragmentMol.delete();
      }
    });
    if (hFragmentSmiles === null)
      throw new Error('Invalid fragment molecule');
    if (!hFragmentSmiles)
      throw new Error('Invalid fragment molecule after adding hydrogens');
    const reactionSmirks = `${hFragmentSmiles}>>[*:1]`;
    const reaction = EngineReaction.create('deprotection', rdModule, reactionSmirks);
    if (!reaction)
      throw new Error('Failed to create reaction from fragment');
    const molList = molColumn.toList();
    // in Crux mode the column on the Crux workers; the molecules they leave, and the column where they fail, here
    const byWorkers = deprotectOnWorkers && reaction.onCrux ? await deprotectByWorkers(molList, reactionSmirks) : null;
    const deprotectedSmiles = molList.map((molSmiles, i) => {
      if (byWorkers?.[i] != null)
        return byWorkers[i];
      try {
        if (molSmiles && !hasNewLines(molSmiles) && molSmiles.length > MAX_SMILES_LENGTH)
          return molSmiles; // do not attempt to parse very long SMILES, will cause MOB.
        return reaction.run((module, rxn) => deprotectOne(module, rxn, molSmiles));
      } catch (_err) {
        return molSmiles;
      }
    });
    reaction.delete();
    const newCol = table.columns.addNewString(table.columns.getUnusedName(`deprotected(${molColumn.name})`));
    newCol.semType = DG.SEMTYPE.MOLECULE;
    newCol.init((i) => deprotectedSmiles[i]);
    return newCol;
  } catch (err: any) {
    grok.shell.error(`Error during deprotection, check console for details`);
    console.error(err);
  }
  return null;
}
