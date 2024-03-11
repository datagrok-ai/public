/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SvgDisplayManager} from './svg-utils/svg-display-manager';

import {EventBus} from '../model/event-bus';
import {PatternEditorDialog} from './pattern-editor';

import $ from 'cash-dom';

export class PatternAppRightSection {
  private svgDisplay: HTMLDivElement;

  constructor(
    private eventBus: EventBus,
  ) {
    this.svgDisplay = SvgDisplayManager.createSvgDiv(eventBus);
  };

  getLayout(): HTMLDivElement {
    const numericLabelTogglesContainer = new NumericLabelVisibilityControls(this.eventBus).getContainer();
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

class NumericLabelVisibilityControls {
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
      const initialValue = this.eventBus.getModificationsWithNumericLabels().includes(nucleotide);
      const toggle = ui.boolInput(
        nucleotide,
        initialValue,
        (value: boolean) => this.handleNumericLabelToggle(nucleotide, value)
      );
      toggles.push(toggle.root);
    });
    return toggles;
  }

  private handleNumericLabelToggle(nucleotide: string, isVisible: boolean): void {
    const labelledNucleotides = this.eventBus.getModificationsWithNumericLabels();
    const hasNumericLabel = labelledNucleotides.includes(nucleotide);
    if (hasNumericLabel === isVisible)
      return;

    const newArray = isVisible ? labelledNucleotides.concat(nucleotide) : labelledNucleotides.filter((n) => n !== nucleotide);
    this.eventBus.updateModificationsWithNumericLabels(newArray);
  }
}
