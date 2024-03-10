/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SvgDisplayManager} from './svg-utils/svg-display-manager';
import _ from 'lodash';

import {EventBus} from '../model/event-bus';
import {BooleanInput} from './types';
import {PatternConfiguration, PhosphorothioateLinkageFlags} from '../model/types';

import $ from 'cash-dom';
import {STRAND, STRANDS} from '../model/const';

export class PatternAppRightSection {
  private svgDisplay: HTMLDivElement;

  constructor(
    private eventBus: EventBus,
  ) {
    this.svgDisplay = SvgDisplayManager.createSvgDiv(eventBus);
  };

  getLayout(): HTMLDivElement {
    const numericLabelTogglesContainer = new NumericLabelToggles(this.eventBus).getContainer();
    const downloadAndEditControls = this.generateDownloadAndEditControls();
    const layout = ui.panel([
      this.svgDisplay,
      numericLabelTogglesContainer,
      downloadAndEditControls,
      // generateStrandSectionDisplays(),
      // ui.h1('Additional modifications'),
      // ui.form([
      //   terminalModificationInputs[STRAND.SENSE][FIVE_PRIME_END],
      //   terminalModificationInputs[STRAND.SENSE][THREE_PRIME_END],
      // ]),
      // asModificationDiv,
      ], {style: {overflowX: 'scroll', padding: '12px 24px'}});
    return layout;
  }

  private generateDownloadAndEditControls(): HTMLDivElement {
    const svgDownloadButton = ui.button('Save PNG', () => this.eventBus.requestSvgSave());

    const editPatternButton = ui.button(
      'Edit pattern',
      () => new PatternEditorDialog(this.eventBus).open()
    );

    return ui.divH([
      svgDownloadButton,
      editPatternButton,
    ], {style: {gap: '12px', marginTop: '12px'}});
  }
}

class NumericLabelToggles {

  private togglesContainer: HTMLDivElement = ui.div([]);

  constructor(
    private eventBus: EventBus,
  ) {
    this.eventBus.uniqueNucleotideBasesChanged$().subscribe(() => {
      this.updateContainer();
    });
  }

  getContainer(): HTMLDivElement {
    return this.togglesContainer;
  }

  private updateContainer(): void {
    $(this.togglesContainer).empty();
    const newToggles = this.createNucleotideToggles();
    this.togglesContainer.append(...newToggles);
  }

  private createNucleotideToggles(): HTMLElement[] {
    const toggles = [] as HTMLElement[];
    const uniqueNucleotideBases = this.eventBus.getUniqueNucleotideBases();
    uniqueNucleotideBases.forEach((nucleotide: string) => {
      const toggle = ui.boolInput(nucleotide, false, (checked: boolean) => {
        grok.shell.info(`Nucleotide ${nucleotide} is ${checked ? 'checked' : 'unchecked'}`);
      });
      toggles.push(toggle.root);
    });
    return toggles;
  }
}

class PatternEditorDialog {
  private initialPatternConfig: PatternConfiguration;

  constructor(
    private eventBus: EventBus,
  ) {
    this.initialPatternConfig = _.cloneDeep(this.eventBus.getPatternConfig());
  }

  open(): void {
    this.createDialog().show();
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

    return allPtoActivationInput;
  }

  private addStyleToPtoInput(allPtoActivationInput: BooleanInput): void {
    const label = allPtoActivationInput.captionLabel;
    ui.tooltip.bind(label, 'Activate all phosphothioate linkages');
    label.classList.add('ui-label-right');
    Object.assign(label.style, {
      textAlign: 'left',
      maxWidth: '100px',
      minWidth: '40px',
      width: 'auto'
    });
  }

  private createDialog(): DG.Dialog {
    const dialog = ui.dialog('Edit pattern')
    .add(ui.divV([
      ui.h1('PTO'),
      ui.divH([
        this.createAllPtoActivationInput().root,
        ...this.createFirstPtoInputs().map((input) => input.root),
      ], {style:{gap:'12px'}})
    ]))
    // .add(ui.divH([
    //   modificationSection[STRAND.SENSE],
    //   modificationSection[STRAND.ANTISENSE],
    // ], {style:{gap:'24px'}}))
    .onOK(() => grok.shell.info('Applied'))
    .onCancel(() => this.resetToInitialState());

    // dialog.onClose.subscribe(() => this.resetToInitialState());

    return dialog;
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
      return firstPtoInput;
    });
  }

  private isFirstPtoActive(strand: STRAND): boolean {
    return this.initialPatternConfig.phosphorothioateLinkageFlags[strand][0];
  }

  private resetToInitialState(): void {
    this.eventBus.setPatternConfig(this.initialPatternConfig);
  }
}
