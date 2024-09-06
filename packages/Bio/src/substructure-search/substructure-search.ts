import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getMonomericMols} from '../calculations/monomerLevelMols';
import {updateDivInnerHTML} from '../utils/ui-utils';
import {delay} from '@datagrok-libraries/utils/src/test';
import {TAGS as bioTAGS, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

export const MONOMER_MOLS_COL = 'monomeric-mols';

export const enum MONOMERIC_COL_TAGS {
  MONOMERIC_MOLS = 'monomeric-mols',
  LAST_INVALIDATED_VERSION = 'last-invalidated-version',
  MONOMERS_DICT = 'monomers-dict'
}

const SUBSTR_HELM_COL_NAME = 'substr_helm';

export class SubstructureSearchDialog {
  units: string;
  separator: string;
  inputsDiv: HTMLDivElement;
  substructureInput: DG.InputBase<string>;
  separatorInput: DG.InputBase<string>;
  editHelmLink: HTMLAnchorElement;
  columnsInput: DG.InputBase<DG.Column | null>;
  grid: DG.Grid;
  col: DG.Column;
  dialog: DG.Dialog;

  constructor(columns: DG.Column<string>[]) {
    this.col = columns[0];
    this.createUI();
  }

  editHelmLinkAction(): void {
    updateDivInnerHTML(this.inputsDiv, this.grid.root);
    ui.tools.waitForElementInDom(this.grid.root).then(() => {
      setTimeout(() => {
        this.grid.cell(SUBSTR_HELM_COL_NAME, 0).element.children[0].dispatchEvent(
          new KeyboardEvent('keydown', {key: 'Enter'})
        );
      }, 100);
    });
  }

  updateInputs(): void {
    const selectedInput = this.units === NOTATION.HELM ? ui.divV([this.columnsInput, this.editHelmLink]) :
      this.units === NOTATION.SEPARATOR ? ui.inputs([this.columnsInput, this.substructureInput, this.separatorInput]) :
        ui.inputs([this.columnsInput, this.substructureInput]);

    updateDivInnerHTML(this.inputsDiv, selectedInput);
  }

  updateNotationDiv(): void {
    this.units = this.col.meta.units!;
    this.separator = this.col.getTag(bioTAGS.separator);
    const notationDiv = this.dialog.root.getElementsByClassName('notation-text')[0];
    if (notationDiv)
      notationDiv.textContent = `Notation: ${this.units}`;
  }

  createUI(): void {
    const dataframe = grok.shell.tv.dataFrame;
    this.columnsInput = ui.input.column('Column', {table: dataframe, value: this.col, onValueChanged: (value) => {
      this.col = value;
      this.updateNotationDiv();
      this.updateInputs();
    }, filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE});

    this.substructureInput = ui.input.string('Substructure', {value: ''});

    this.editHelmLink = ui.link('Edit helm', () => this.editHelmLinkAction(), undefined, {style: {position: 'relative', left: '95px'}});

    const df = DG.DataFrame.create(1);
    df.columns.addNewString(SUBSTR_HELM_COL_NAME).init((_i) => '');
    df.col(SUBSTR_HELM_COL_NAME)!.semType = this.col.semType;
    df.col(SUBSTR_HELM_COL_NAME)!.meta.units = NOTATION.HELM;
    this.grid = df.plot.grid();
    this.separatorInput = ui.input.string('Separator', {value: this.separator});

    this.inputsDiv = ui.div();
    this.units = this.col.meta.units!;
    this.separator = this.col.getTag(bioTAGS.separator);
    this.updateInputs();

    this.dialog = ui.dialog('Substructure Search')
      .add(ui.divV([
        ui.divText(`Notation: ${this.units}`, 'notation-text'),
        this.inputsDiv,
      ]))
      .onOK(async () => {
        let substructure = this.units === NOTATION.HELM ? df.get(SUBSTR_HELM_COL_NAME, 0) : this.substructureInput.value;
        if (this.units === NOTATION.SEPARATOR && this.separatorInput.value !== this.separator && this.separatorInput.value !== '')
          substructure = substructure.replaceAll(this.separatorInput.value, this.separator);
        let matches: DG.BitSet;
        if (this.units === NOTATION.HELM)
          matches = await helmSubstructureSearch(substructure, this.col);
        else
          matches = linearSubstructureSearch(substructure, this.col);
        this.col.dataFrame.filter.and(matches);
      })
      .show();
  }
}

export function linearSubstructureSearch(substructure: string, col: DG.Column<string>, separator?: string): DG.BitSet {
  const re = separator ? prepareSubstructureRegex(substructure, separator) : substructure;
  const resultArray = DG.BitSet.create(col.length);
  for (let i = 0; i < col.length; i++) {
    const macromolecule: string = col.get(i)!;
    if (macromolecule.match(re) || macromolecule === substructure)
      resultArray.set(i, true, false);
  }
  return resultArray;
}

function prepareSubstructureRegex(substructure: string, separator: string) {
  const char = `${separator}`.replace(/[\-\[\]\/\{\}\(\)\*\+\?\.\\\^\$\|]/g, '\\$&');
  const startsWithSep = substructure.charAt(0) === separator;
  const endsWithSep = substructure.charAt(substructure.length - 1) === separator;
  const substrWithoutSep = substructure.replace(new RegExp(`^${char}|${char}$`, 'g'), '');
  const re = startsWithSep ?
    endsWithSep ? `${char}${substrWithoutSep}${char}` :
      `${char}${substrWithoutSep}${char}|${char}${substrWithoutSep}$` :
    endsWithSep ? `^${substrWithoutSep}${char}|${char}${substrWithoutSep}${char}` :
      `^${substrWithoutSep}${char}|${char}${substrWithoutSep}${char}|${char}${substrWithoutSep}$`;
  return re;
}

export async function helmSubstructureSearch(substructure: string, col: DG.Column<string>): Promise<DG.BitSet> {
  if (col.version !== col.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION])
    await invalidateMols(col, true);
  const substructureCol: DG.Column<string> = DG.Column.string('helm', 1).init((_i) => substructure);
  substructureCol.semType = DG.SEMTYPE.MACROMOLECULE;
  substructureCol.meta.units = NOTATION.HELM;
  const substructureMolsCol =
    await getMonomericMols(substructureCol, true, col.temp[MONOMERIC_COL_TAGS.MONOMERS_DICT]);
  const matchesCol = await grok.functions.call('Chem:searchSubstructure', {
    molStringsColumn: col.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS],
    molString: substructureMolsCol.get(0),
    molBlockFailover: '',
  });
  return matchesCol.get(0);
}

export async function invalidateMols(col: DG.Column<string>, pattern: boolean) {
  const progressBar = DG.TaskBarProgressIndicator.create(`Invalidating molfiles for ${col.name}`);
  try {
    await delay(10);
    const monomersDict = new Map();
    const monomericMolsCol = await getMonomericMols(col, pattern, monomersDict);
    col.temp[MONOMERIC_COL_TAGS.MONOMERIC_MOLS] = monomericMolsCol;
    col.temp[MONOMERIC_COL_TAGS.MONOMERS_DICT] = monomersDict;
    col.temp[MONOMERIC_COL_TAGS.LAST_INVALIDATED_VERSION] = col.version;
  } finally {
    progressBar.close();
  }
}
