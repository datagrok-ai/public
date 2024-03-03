/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import { STRAND, STRAND_LABEL, STRANDS, OTHER_USERS } from '../model/const';
import { PatternConfiguration, StrandType } from '../model/types';

import {StringInput, NumberInput} from './types';
import {NucleotidePatternSVGRenderer} from '../view/svg-utils/svg-renderer';

import {EventBus} from '../model/event-bus';
import {PatternAppDataManager} from '../model/external-data-manager';
import * as rxjs from 'rxjs';
import $ from 'cash-dom';
import {PatternConfigurationManager} from '../model/pattern-config-manager';

export class PatternAppRightSection {
  constructor(
    private eventBus: EventBus,
    private dataManager: PatternAppDataManager,
  ) { };

  getLayout(): HTMLDivElement {
    const svgDisplay = SvgDisplayManager.createSvgDiv(this.eventBus);
    const layout = ui.panel([
      svgDisplay,
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

class SvgDisplayManager {
  private svgDisplayDiv = ui.div([]);
  private configManager: PatternConfigurationManager;

  private constructor(
    eventBus: EventBus,
  ) {
    this.configManager = new PatternConfigurationManager(eventBus);
  }

  static createSvgDiv(eventBus: EventBus): HTMLElement {
    const manager = new SvgDisplayManager(eventBus);
    manager.updateSvgContainer();
    return manager.svgDisplayDiv;
  }

  private updateSvgContainer(): void {
    $(this.svgDisplayDiv).empty();
    const patternConfig = this.configManager.getConfig();
    const image = this.createSvg(patternConfig);
    this.svgDisplayDiv.append(image);
  }

  private createSvg(patternConfig: PatternConfiguration) {
    const renderer = new NucleotidePatternSVGRenderer(patternConfig);
    const svg = renderer.renderPattern();
    return svg;
  }
}
