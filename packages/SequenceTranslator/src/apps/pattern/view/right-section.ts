/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SvgDisplayManager} from './svg-utils/svg-display-manager';

import {EventBus} from '../model/event-bus';
import {BooleanInput} from './types';
import {PatternConfigurationManager} from '../model/pattern-config-manager';
import {PatternConfiguration} from '../model/types';

import $ from 'cash-dom';
import {STRANDS} from '../model/const';

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
  private patternConfig: PatternConfiguration;

  constructor(
    private eventBus: EventBus,
  ) {
    this.patternConfig = PatternConfigurationManager.getConfig(this.eventBus);
  }

  open(): void {
    this.createDialog().show();
  }

  private createAllPtoActivationInput(): BooleanInput {
    const flags = this.patternConfig.phosphorothioateLinkageFlags;
    const totalNumberOfPTOFlags = STRANDS.map((strand) => flags[strand].length).reduce((a, b) => a + b, 0);
    const totalNumberOfNucleotides = STRANDS.map((strand) => this.patternConfig.nucleotideSequences[strand].length).reduce((a, b) => a + b, 0);

    // WARNING: +2 because there are +1 more PTO flags in each strand than there are nucleotides
    const allPTOLinkagesSet = totalNumberOfPTOFlags === totalNumberOfNucleotides + 2;

    const allPtoActivationInput = ui.boolInput('All PTO',
      allPTOLinkagesSet,
      (value: boolean) => this.eventBus.setAllPTOLinkages(value)
    );

    ui.tooltip.bind(allPtoActivationInput.root, 'Activate all phosphothioate linkages');

    const label = allPtoActivationInput.captionLabel;

    label.classList.add('ui-label-right');
    Object.assign(label.style, {
      textAlign: 'left',
      maxWidth: '100px',
      minWidth: '40px',
      width: 'auto'
    });
    return allPtoActivationInput;
  }

  private createDialog(): DG.Dialog {
    return ui.dialog('Edit pattern')
    .add(ui.divV([
      ui.h1('PTO'),
      ui.divH([
        this.createAllPtoActivationInput().root,
        // firstPto[STRAND.SENSE].root,
        // firstPto[STRAND.ANTISENSE].root,
      ], {style:{gap:'12px'}})
    ]))
    // .add(ui.divH([
    //   modificationSection[STRAND.SENSE],
    //   modificationSection[STRAND.ANTISENSE],
    // ], {style:{gap:'24px'}}))
    .onOK(()=>{grok.shell.info('Saved')})
  }
}
