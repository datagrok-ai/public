/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SvgDisplayManager} from './svg-utils/svg-display-manager';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import {PatternConfigurationManager} from '../model/pattern-config-manager';

export class PatternAppRightSection {
  private svgDisplay: HTMLDivElement;

  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
    patternConfiguration: PatternConfigurationManager,
  ) {
    this.svgDisplay = SvgDisplayManager.createSvgDiv(eventBus, patternConfiguration);
  };

  getLayout(): HTMLDivElement {
    const layout = ui.panel([
      this.svgDisplay,
      // numericLabelTogglesContainer,
      // generateDownloadAndEditControls(),
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

  // private createNumericLabelTogglesContainer(): HTMLElement {
  //   const defaultNucleotideBase = this.patternConfiguration.getDefaultNucleotideBase();
  //   const numericLabelTogglesContainer = ui.divH([
  //     ui.boolInput(defaultNucleotideBase, true, (value: boolean) => {
  //       updateListOfModificationsWithNumericLabels(value);
  //       refreshSvgDisplay();
  //       refreshOutputExamples();
  //     }).root,
  //   ]);
  //   return numericLabelTogglesContainer;
  // }
}

