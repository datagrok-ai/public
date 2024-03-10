/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import _ from 'lodash';

import {EventBus} from '../model/event-bus';
import {BooleanInput} from './types';
import {PatternConfiguration, PhosphorothioateLinkageFlags} from '../model/types';
import {STRAND, STRANDS} from '../model/const';

export class PatternEditorDialog {
  private initialPatternConfig: PatternConfiguration;

  constructor(
    private eventBus: EventBus,
  ) {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
  }

  open(): void {
    this.createDialog().show();
  }

  private createDialog(): DG.Dialog {
    const header = new HeaderControls(this.eventBus, this.initialPatternConfig).getPhosphorothioateLinkageControls();
    const dialog = ui.dialog('Edit pattern')
    .add(header)
    // .add(ui.divH([
    //   modificationSection[STRAND.SENSE],
    //   modificationSection[STRAND.ANTISENSE],
    // ], {style:{gap:'24px'}}))
    .onOK(() => grok.shell.info('Applied'))
    .onCancel(() => this.resetToInitialState());

    // dialog.onClose.subscribe(() => this.resetToInitialState());

    return dialog;
  }

  private resetToInitialState(): void {
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}

class HeaderControls {
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
      ], {style:{gap:'12px'}})
    ])

    return container;
  }

  private areAllPtoLinkagesSet(flags: PhosphorothioateLinkageFlags): boolean {
    const totalNumberOfPTOFlags = STRANDS.map(
      (strand) => flags[strand].filter((flag) => flag).length
    )
      .reduce((a, b) => a + b, 0);
    const totalNumberOfNucleotides = STRANDS.map((strand) => this.initialPatternConfig.nucleotideSequences[strand].length).reduce((a, b) => a + b, 0);

    // WARNING: +2 because there are +1 more PTO flags in each strand than there are nucleotides
    return totalNumberOfPTOFlags === totalNumberOfNucleotides + 2;
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
    ui.tooltip.bind(allPtoActivationInput.captionLabel, 'Include all phosphothioates in the pattern');

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

  private createFirstPtoInputs(): BooleanInput[] {
    return STRANDS.map((strand) => {
      const initialValue = this.isFirstPtoActive(strand);
      const firstPtoInput = ui.boolInput(`First ${strand} PTO`, initialValue, );

      firstPtoInput.onInput(() => {
        const value = firstPtoInput.value!;
        this.eventBus.setFirstPhosphorothioateLinkageFlag(strand, value);
      });

      this.eventBus.phosphorothioateLingeFlagsChanged$.subscribe((flags) => {
        const newValue = flags[strand][0];
        firstPtoInput.value = newValue;
      });
      
      this.addStyleToPtoInput(firstPtoInput);
      ui.tooltip.bind(firstPtoInput.captionLabel, `Include first phosphothioate in ${strand}`);
      return firstPtoInput;
    });
  }

  private isFirstPtoActive(strand: STRAND): boolean {
    return this.initialPatternConfig.phosphorothioateLinkageFlags[strand][0];
  }
}

class StrandControls {
  constructor(
    private eventBus: EventBus,
    private initialPatternConfig: PatternConfiguration,
  ) { }
}
