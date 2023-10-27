/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';
import '../style/translator-app.css';
import $ from 'cash-dom';

import {highlightInvalidSubsequence} from '../utils/colored-input/input-painters';
import {ColoredTextInput} from '../utils/colored-input/colored-text-input';
import {SequenceToMolfileConverter} from '../../model/sequence-to-structure-utils/sequence-to-molfile';
import {getTranslatedSequences} from '../../model/format-translation/conversion-utils';
import {MoleculeImage} from '../utils/molecule-img';
import {download} from '../../model/helpers';
import {SEQUENCE_COPIED_MSG, SEQ_TOOLTIP_MSG} from '../const/oligo-translator';
import {DEFAULT_AXOLABS_INPUT} from '../const/ui';
import {FormatDetector} from '../../model/parsing-validation/format-detector';
import {SequenceValidator} from '../../model/parsing-validation/sequence-validator';
import {FormatConverter} from '../../model/format-translation/format-converter';
import {codesToHelmDictionary} from '../../model/data-loading-utils/json-loader';
import {DEFAULT_FORMATS} from '../../model/const';

export class TranslatorLayoutHandler {
  constructor() {
    const INPUT_FORMATS = Object.keys(codesToHelmDictionary).concat(DEFAULT_FORMATS.HELM);
    this.moleculeImgDiv = ui.div([]);
    this.moleculeImgDiv.className = 'mol-host';
    this.moleculeImgDiv.style.border = '1px solid var(--grey-2)';
    this.moleculeImgDiv.style.borderRadius = '1px';
    this.moleculeImgDiv.style.marginTop = '12px';

    this.outputTableDiv = ui.div([]);
    this.formatChoiceInput = ui.choiceInput('', DEFAULT_FORMATS.HELM, INPUT_FORMATS, async () => {
      this.format = this.formatChoiceInput.value;
      this.updateTable();
      await this.updateMolImg();
    });
    this.sequenceInputBase = ui.textInput('', DEFAULT_AXOLABS_INPUT,
      () => { this.onInput.next(); });

    this.init();

    DG.debounce<string>(this.onInput, 300).subscribe(async () => {
      this.init();
      this.formatChoiceInput.value = this.format;
      this.updateTable();
      await this.updateMolImg();
    });
  }

  // todo: reduce # of state vars by further refactoring legacy code
  private onInput = new rxjs.Subject<string>();
  private moleculeImgDiv: HTMLDivElement;
  private outputTableDiv: HTMLDivElement;
  private formatChoiceInput: DG.InputBase;
  private sequenceInputBase: DG.InputBase;
  private molfile: string;
  public sequence: string;
  private format: string | null;

  async getHtmlElement(): Promise<HTMLDivElement> {
    const sequenceColoredInput = new ColoredTextInput(this.sequenceInputBase, highlightInvalidSubsequence);

    const downloadMolfileButton = ui.button(
      'Get SDF',
      () => { this.saveMolfile(); },
      'Save structure as SDF');

    const copySmilesButton = ui.button(
      'Copy SMILES',
      () => { this.copySmiles(); },
      'Copy SMILES for the sequence');

    const formatChoiceInput = ui.div([this.formatChoiceInput]);

    const inputTableRow = {
      format: formatChoiceInput,
      textInput: sequenceColoredInput.root,
    };
    const upperBlock = ui.table(
      [inputTableRow], (item) => [item.format, item.textInput]
    );
    upperBlock.classList.add('st-translator-input-table');

    const outputTable = ui.block([
      this.outputTableDiv,
      downloadMolfileButton,
      copySmilesButton,
    ]);

    const mainTabBody = ui.box(
      ui.panel([
        upperBlock,
        outputTable,
        ui.block([ui.box(this.moleculeImgDiv)])
      ]),
    );

    this.formatChoiceInput.value = this.format;
    this.updateTable();
    await this.updateMolImg();
    return mainTabBody;
  }

  private saveMolfile(): void {
    const result = (new SequenceToMolfileConverter(this.sequence, false,
      this.formatChoiceInput.value!)).convert() + '\n$$$$';
    download(this.sequence + '.sdf', encodeURIComponent(result));
  }

  private copySmiles(): void {
    const smiles = DG.chem.convert(this.molfile, DG.chem.Notation.MolBlock, DG.chem.Notation.Smiles);
    navigator.clipboard.writeText(smiles).then(
      () => grok.shell.info(SEQUENCE_COPIED_MSG)
    );
  }

  private updateTable(): void {
    this.outputTableDiv.innerHTML = '';
    // todo: does not detect correctly (U-A)(U-A)
    const indexOfInvalidChar = (!this.format) ? 0 : (new SequenceValidator(this.sequence)).getInvalidCodeIndex(this.format!);
    const translatedSequences = getTranslatedSequences(this.sequence, indexOfInvalidChar, this.format!);
    const tableRows = [];

    for (const key of Object.keys(translatedSequences)) {
      const sequence = ('indexOfFirstInvalidChar' in translatedSequences) ?
        ui.divH([]) :
        ui.link(
          translatedSequences[key],
          () => navigator.clipboard.writeText(translatedSequences[key])
            .then(() => grok.shell.info(SEQUENCE_COPIED_MSG)),
          SEQ_TOOLTIP_MSG, ''
        );
      tableRows.push({
        format: key,
        sequence: sequence,
      });
    }
    const outputTable = ui.table(tableRows, (item) => [item.format, item.sequence], ['FORMAT', 'SEQUENCE']);

    this.outputTableDiv.append(outputTable);
    this.outputTableDiv.classList.add('st-translator-output-table');
  }

  private async updateMolImg(): Promise<void> {
    const canvasWidth = 500;
    const canvasHeight = 170;
    const molImgObj = new MoleculeImage(this.molfile);
    await molImgObj.drawMolecule(this.moleculeImgDiv, canvasWidth, canvasHeight);
    // should the canvas be returned from the above function?
  }

  // todo: sort mehtods
  private init(): void {
    this.sequence = this.getFormattedSequence();

    this.format = (new FormatDetector(this.sequence)).getFormat();

    // warning: getMolfile relies on 'this.format', so the order is important
    this.molfile = this.getMolfile();
  }

  private getFormattedSequence(): string {
    return this.sequenceInputBase.value.replace(/\s/g, '');
  }

  private getMolfile(): string {
    if (!this.format)
      return '';
    if (this.format === DEFAULT_FORMATS.HELM) {
      const axolabs = (new FormatConverter(this.sequence, this.format).convertTo(DEFAULT_FORMATS.AXOLABS));
      return (new SequenceToMolfileConverter(axolabs, false, DEFAULT_FORMATS.AXOLABS).convert());
    }
    const molfile = (new SequenceToMolfileConverter(this.sequence, false, this.format)).convert();
    return molfile;
  }
}
