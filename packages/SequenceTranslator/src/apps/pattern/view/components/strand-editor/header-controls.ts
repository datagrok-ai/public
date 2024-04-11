import * as ui from 'datagrok-api/ui';

import {STRAND, STRANDS} from '../../../model/const';
import {EventBus} from '../../../model/event-bus';
import {PatternConfiguration, PhosphorothioateLinkageFlags} from '../../../model/types';
import {BooleanInput} from '../../types';

export class HeaderControls {
  constructor(
    private eventBus: EventBus,
    private initialPatternConfig: PatternConfiguration
  ) { }

  getPhosphorothioateLinkageControls(): HTMLDivElement {
    const container = ui.divV([
      ui.h1('PTO'),
      ui.divH([
        this.createAllPtoActivationInput().root,
        ...this.createFirstPtoInputs().map((input) => input.root),
      ], {style: {gap: '12px'}})
    ]);

    return container;
  }

  private areAllPtoLinkagesSet(flags: PhosphorothioateLinkageFlags): boolean {
    const totalNumberOfPTOFlags = STRANDS.map(
      (strand) => flags[strand].filter((flag) => flag).length
    )
      .reduce((a, b) => a + b, 0);
    const totalNumberOfNucleotides = STRANDS.map(
      (strand) => this.initialPatternConfig.nucleotideSequences[strand].length
    ).reduce((a, b) => a + b, 0);

    // There are +1 more PTO flags in each strand than there are nucleotides
    const addendum = STRANDS.filter((strand) => flags[strand].length).length;
    return totalNumberOfPTOFlags === totalNumberOfNucleotides + addendum;
  }

  private createAllPtoActivationInput(): BooleanInput {
    const flags = this.initialPatternConfig.phosphorothioateLinkageFlags;
    const initialValue = this.areAllPtoLinkagesSet(flags);
    const allPtoActivationInput = ui.boolInput('All PTO', initialValue);

    allPtoActivationInput.onInput(() => {
      const value = allPtoActivationInput.value!;
      this.eventBus.setAllPTOLinkages(value);
    });

    this.eventBus.phosphorothioateLingeFlagsChanged$.subscribe(() => {
      const flags = this.eventBus.getPatternConfig().phosphorothioateLinkageFlags;
      const newValue = this.areAllPtoLinkagesSet(flags);
      allPtoActivationInput.value = newValue;
    });

    this.addStyleToPtoInput(allPtoActivationInput);
    ui.tooltip.bind(allPtoActivationInput.captionLabel, 'Activate all phosphothioates');

    return allPtoActivationInput;
  }

  private addStyleToPtoInput(allPtoActivationInput: BooleanInput): void {
    const label = allPtoActivationInput.captionLabel;
    label.classList.add('ui-label-right');
    Object.assign(label.style, {
      textAlign: 'left',
      maxWidth: '100px',
      minWidth: '40px',
      width: 'auto'
    });
  }

  private createFirstPtoInputs(): BooleanInput [] {
    return STRANDS.map((strand) => {
      if (!this.eventBus.isAntisenseStrandActive() && strand === STRAND.ANTISENSE)
        return;
      const initialValue = this.isFirstPtoActive(strand);
      const firstPtoInput = ui.boolInput(`First ${strand} PTO`, initialValue);

      firstPtoInput.onInput(() => {
        const value = firstPtoInput.value!;
        this.eventBus.setPhosphorothioateLinkageFlag(strand, 0, value);
      });

      this.eventBus.phosphorothioateLingeFlagsChanged$.subscribe((flags) => {
        const newValue = flags[strand][0];
        firstPtoInput.value = newValue;
      });

      this.addStyleToPtoInput(firstPtoInput);
      ui.tooltip.bind(firstPtoInput.captionLabel, `Activate first phosphothioate in ${strand}`);
      return firstPtoInput;
    }).filter((input) => input !== undefined) as BooleanInput[];
  }

  private isFirstPtoActive(strand: STRAND): boolean {
    return this.initialPatternConfig.phosphorothioateLinkageFlags[strand][0];
  }
}

