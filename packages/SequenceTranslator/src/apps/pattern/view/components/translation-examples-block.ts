import * as ui from 'datagrok-api/ui';

import '../style.css';

import {StringInput} from '../types';

import $ from 'cash-dom';
import {NUCLEOTIDES} from '../../../common/model/const';
import {STRAND, STRANDS, STRAND_LABEL} from '../../model/const';
import {DataManager} from '../../model/data-manager';
import {EventBus} from '../../model/event-bus';
import {applyPatternToRawSequence} from '../../model/translator';


export class TranslationExamplesBlock {
  constructor(
    private eventBus: EventBus,
    private dataManager: DataManager
  ) { }

  createContainer(): HTMLDivElement {
    return ui.div([
      ui.h1('Translation example'),
      this.createTranslationExamples(),
    ], {style: {paddingTop: '20px'}});
  }

  private createTranslationExamples(): HTMLDivElement {
    const strandExamples = ui.divH([
      ...this.getExampleElements(),
    ], 'ui-form');

    this.eventBus.strandsAndLinkagesUpdated$.subscribe(() => {
      $(strandExamples).empty();
      $(strandExamples).append(this.getExampleElements());
    });

    return strandExamples;
  }

  private getExampleElements(): HTMLDivElement[] {
    return STRANDS.map((strand) => this.createStrandExample(strand));
  }

  private createStrandExample(strand: STRAND): HTMLDivElement {
    if (!this.eventBus.isAntisenseStrandActive() && strand === STRAND.ANTISENSE)
      return ui.div([]);

    const inputExample = this.createInputExample(strand);
    const rawSequence = inputExample.value!;
    const outputExample = this.createOutputExample(strand, rawSequence);
    const strandExample = ui.block50([
      ui.h2(STRAND_LABEL[strand]),
      inputExample.root,
      outputExample.root,
    ], {style: {paddingRight: '20px'}});
    return strandExample;
  }

  private createInputExample(strand: STRAND): StringInput {
    const input = this.createTextInputForExamples();

    const exampleRawNucleotides = this.generateExampleSequence(strand);
    input.value = exampleRawNucleotides;

    input.setTooltip(`Example raw nucleotides input for ${STRAND_LABEL[strand]}`);

    return input;
  }

  private generateExampleSequence(strand: STRAND): string {
    const sourceSequence = this.eventBus.getNucleotideSequences()[strand];
    const exampleSequence = sourceSequence.map((_, index) => {
      return NUCLEOTIDES[index % NUCLEOTIDES.length];
    }).join('');

    return exampleSequence;
  }

  private createOutputExample(strand: STRAND, rawNucleotideSequence: string): StringInput {
    const output = this.createTextInputForExamples();

    const modifications = this.eventBus.getNucleotideSequences()[strand];
    const ptoFlags = this.eventBus.getPhosphorothioateLinkageFlags()[strand];
    output.value = applyPatternToRawSequence(rawNucleotideSequence, modifications, ptoFlags);

    output.setTooltip(`Pattern applied to the example input for ${STRAND_LABEL[strand]}`);

    return output;
  }

  private createTextInputForExamples(): StringInput {
    const input = ui.textInput('', '');
    this.applyStylingToInput(input);

    return input;
  }

  private applyStylingToInput(input: StringInput): void {
    const textarea = input.root.getElementsByTagName('textarea')[0];
    textarea.setAttribute('readonly', 'true');
    $(textarea).css('resize', 'none');
    $(input.root).css('opacity', '75%');
  }
}
