import * as ui from 'datagrok-api/ui';

import {SubscriptionManager} from '../../model/subscription-manager';

import '../style.css';

import {StringInput} from '../types';

import $ from 'cash-dom';
import {DataManager} from '../../model/data-manager';
import {EventBus} from '../../model/event-bus';
import {STRAND, STRANDS, STRAND_LABEL} from '../../model/const';
import {NUCLEOTIDES} from '../../../common/model/const';


export class TranslationExamplesBlock {
  private subscriptions = new SubscriptionManager();

  constructor(
    private eventBus: EventBus,
    private dataManager: DataManager
  ) { }

  createContainer(): HTMLDivElement {
    return ui.div([
      ui.h1('Translation examples'),
      this.createStrandExamples(),
    ], {style: {paddingTop: '20px'}});
  }

  private createStrandExamples(): HTMLDivElement {
    const strandExamples = ui.divH([
      ...STRANDS.map((strand) => this.createStrandExample(strand)),
    ], 'ui-form');

    return strandExamples;
  }

  private createStrandExample(strand: STRAND): HTMLDivElement {
    const inputExample = this.createInputExample(strand);
    const outputExample = this.createOutputExample(strand);
    const strandExample = ui.block([
      ui.h2(STRAND_LABEL[strand]),
      inputExample.root,
      outputExample.root,
    ], {style: {paddingRight: '20px'}});
    return strandExample;
  }

  private createInputExample(strand: STRAND): StringInput {
    const exampleRawNucleotides = this.generateExampleSequence(strand);
    const input = ui.textInput('', exampleRawNucleotides);

    const textarea = input.root.getElementsByTagName('textarea')[0];
    textarea.setAttribute('readonly', 'true');
    $(textarea).css('resize', 'none');
    $(input.root).css('opacity', '75%');

    input.setTooltip(`Example input for ${STRAND_LABEL[strand]}`);

    return input;
  }

  private generateExampleSequence(strand: STRAND): string {
    const sourceSequence = this.eventBus.getNucleotideSequences()[strand];
    const exampleSequence = sourceSequence.map((_, index) => {
      return NUCLEOTIDES[index % NUCLEOTIDES.length];
    }).join('');

    return exampleSequence;
  }

  private createOutputExample(strand: STRAND): StringInput {
    return this.createInputExample(strand);
  }
}
