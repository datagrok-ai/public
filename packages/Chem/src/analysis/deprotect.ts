/* eslint-disable max-len */
import {MolList, RDModule, RDMol, RDReactionResult} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Observable, Subject} from 'rxjs';
import {hasNewLines} from '../utils/chem-common';

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

const defaultFramgment = 'O=C(N[*:1])OCC1c2ccccc2-c2ccccc21';
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
    this.fragmentInput = ui.input.molecule('Fragment', {value: radicalToRgroupForm(funcCall.inputs['fragment'] ?? defaultFramgment), nullable: false, onValueChanged: () => {
      funcCall.inputs['fragment'] = correctRGroups(grok.chem.convert(this.fragmentInput.value ?? '', grok.chem.Notation.Unknown, grok.chem.Notation.Smiles));
      this.inputChangedSubject.next();
    }, tooltipText: 'Smiles of the protecting group fragment to be removed. R1 group should be used to mark the attachment point.'});
    root.appendChild(this.tableInput.root);
    root.appendChild(this.fragmentInput.root);
    this.updateColumnInput();
    this.fragmentInput.addValidator(() => this.validate());
    this.tableInput.addValidator(() => this.validate());
    // triger changes to init the funccall
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
    if (correctedFragment.indexOf('[*:1]') == -1)
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

export function addDeprotectedColumn(table: DG.DataFrame, molColumn: DG.Column, fragment: string, rdModule: RDModule): DG.Column | null {
  try {
    const fragmentMol = rdModule.get_mol(fragment);
    if (!fragmentMol)
      throw new Error('Invalid fragment molecule');
    // we need to add explicit hydrogens to the fragment to make sure the reaction works correctly
    fragmentMol.add_hs_in_place();
    const hFragmentSmiles = fragmentMol.get_smiles();
    fragmentMol.delete();
    if (!hFragmentSmiles)
      throw new Error('Invalid fragment molecule after adding hydrogens');
    const reactionSmirks = `${hFragmentSmiles}>>[*:1]`;
    const reaction = rdModule.get_rxn(reactionSmirks);
    if (!reaction)
      throw new Error('Failed to create reaction from fragment');
    const molList = molColumn.toList();
    const deprotectedSmiles = molList.map((molSmiles) => {
      let mols: MolList | null = null;
      let mol: RDMol | null = null;
      let rctns: RDReactionResult | null = null;
      try {
        if (molSmiles && !hasNewLines(molSmiles) && molSmiles.length > 5000)
          return molSmiles; // do not attempt to parse very long SMILES, will cause MOB.
        mol = rdModule.get_mol(molSmiles);
        if (!mol)
          return molSmiles;
        mol.add_hs_in_place();
        mols = new rdModule.MolList();
        mols.append(mol);
        rctns = reaction.run_reactants(mols, 1);
        const size = rctns.size();
        if (size === 0)
          return molSmiles;
        const firstProduct = rctns.get(0).next();
        if (!firstProduct)
          return molSmiles;
        firstProduct.remove_hs_in_place();
        const deprotectedBlock = firstProduct.get_molblock();
        return deprotectedBlock ?? molSmiles;
      } catch (_err) {
        return molSmiles;
      } finally {
        mols?.delete();
        mol?.delete();
        rctns?.delete();
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
