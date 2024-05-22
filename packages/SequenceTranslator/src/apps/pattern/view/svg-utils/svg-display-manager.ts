/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PatternConfiguration, StrandType} from '../../model/types';
import {EventBus} from '../../model/event-bus';
import {NucleotidePatternSVGRenderer} from './svg-renderer';
//@ts-ignore
import * as svgExport from 'save-svg-as-png';

export class SvgDisplayManager {
  private svgDisplayDiv = ui.div([]);
  private svgElement: SVGElement;

  private constructor(
    private eventBus: EventBus
  ) {
    eventBus.patternStateChanged$.subscribe(() => this.updateSvgContainer());
    eventBus.svgSaveRequested$.subscribe(() => this.saveSvgAsPng());
  }

  static createSvgDiv(eventBus: EventBus): HTMLDivElement {
    const displayManager = new SvgDisplayManager(eventBus);
    return displayManager.svgDisplayDiv;
  }

  private updateSvgContainer(): void {
    $(this.svgDisplayDiv).empty();
    const patternConfig = this.eventBus.getPatternConfig();
    this.svgElement = this.createSvg(patternConfig);
    this.svgDisplayDiv.append(this.svgElement);
  }

  private createSvg(patternConfig: PatternConfiguration) {
    const renderer = new NucleotidePatternSVGRenderer(patternConfig, this.eventBus);
    const svg = renderer.renderPattern();
    return svg;
  }

  private saveSvgAsPng(): void {
    const patternName = this.eventBus.getPatternName();
    svgExport.saveSvgAsPng(this.svgElement, patternName, {backgroundColor: 'white'});
  }
}
