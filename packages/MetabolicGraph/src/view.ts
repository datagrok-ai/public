import * as DG from 'datagrok-api/dg';

import type {BuilderType, BuilderConstructor} from '../escher_src/src/Builder';
import type {MapData, CobraModelData} from '../escher_src/src/ts/types';
import {sampleReactions, saveStateDialog, loadAnalisisDialog, runFBADialog} from './utils';
import {clearTimeCourseSlider} from './timeCourse';

/** Handed to the registered AI view functions (package.ts / ai-tools.ts). */
export interface MetabolicAIContext {
  builder: BuilderType;
}

export class MetabolicGraphView extends DG.ViewBase {
  builder: BuilderType | null = null;

  constructor() {
    super();
    this.name = 'Metabolic Graph';
    this.root.classList.add('d4-escher-container');
    this.aiDescription = 'Metabolic Graph — an Escher-based metabolic network map backed by a COBRA model, ' +
      'for flux analysis (FBA, flux sampling, time-course simulation). Act on it through the view functions ' +
      '(search list_view_functions with "metabolic"): getMetabolicState (call first to see what is loaded), ' +
      'findMetabolicEntities / getMetabolicEntityDetails to explore reactions, metabolites and genes, ' +
      'setReactionBounds / runFba / sampleReactionFluxes / runTimeCourse to run analyses (results color the map), ' +
      'setReactionFluxes to color the map with custom data, showMetabolicEntity / highlightShortestPath to point ' +
      'things out to the user, and listMetabolicAnalyses / saveMetabolicAnalysis / loadMetabolicAnalysis for ' +
      'saved analyses.';
  }

  /** Collected by the AI assistant through `view.getFunctions()` (the Dart JsViewHost forwards here). */
  getFunctions(): DG.Func[] {
    return DG.Func.find({package: 'MetabolicGraph', tags: ['metabolicViewFunction']});
  }

  /** Facade the registered AI view functions use to act on this instance. */
  aiContext(): MetabolicAIContext | null {
    return this.builder ? {builder: this.builder} : null;
  }
}

/** Builds the Escher builder inside [view] with the standard Datagrok wiring (sampling,
 * save/load, FBA dialogs, time-course cleanup) and registers it on the view. */
export function createEscherBuilder(view: MetabolicGraphView, mapData: MapData,
  modelData: CobraModelData | null): BuilderType {
  const Builder = window.escher.Builder as BuilderConstructor;
  const b: BuilderType = new Builder(mapData, modelData, null, window.escher.libs.d3_select(view.root),
    {scroll_behavior: 'zoom', fill_screen: false,
      never_ask_before_quit: true,
      samplingFunction: (mp: CobraModelData) => sampleReactions(mp, b),
      saveAction: () => saveStateDialog(b, undefined),
      loadAction: () => loadAnalisisDialog(b),
      runFBA: async () => { runFBADialog(b); },
      pathFindingDisabled: true,
    });
  // remove the time-course animation slider whenever the map/model is replaced or cleared
  b.callback_manager?.set('load_map.timecourse', () => clearTimeCourseSlider());
  b.callback_manager?.set('load_model.timecourse', () => clearTimeCourseSlider());
  b.callback_manager?.set('clear_map.timecourse', () => clearTimeCourseSlider());
  view.builder = b;
  return b;
}
