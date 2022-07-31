import {HitTriageBaseView} from "./hit-triage-base-view";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {HitTriageApp} from "../hit-triage-app";


/**
 * Enrichment of the molecular dataset.
 **/
export class EnrichView extends HitTriageBaseView {

  grid: DG.Grid;

  constructor(app: HitTriageApp) {
    super(app);

    this.grid = this.template.sourceDataFrame!.plot.grid()
    const content = ui.divV([
      ui.h1('This is where we enrich our data'),
      this.grid
    ])

    this.root.appendChild(content);
  }

  onActivated() {
    super.onActivated();
    this.process().then(() => {
    });
  }

  async process(): Promise<any> {
    this.template.enrichedDataFrame = this.template.sourceDataFrame;

    await this.template.enrichedDataFrame!.columns.addNewCalculated("length", "length(${smiles})");

    // descriptors
    let descriptors = ["HeavyAtomCount", "NHOHCount"];
    await grok.chem.descriptors(this.template.enrichedDataFrame!, this.template.sourceMoleculeColumn, descriptors);

    // another way to invoke calculations
    //await grok.functions.call('ChemDescriptors', {table: session.enrichedDataFrame, smiles: session.sourceMoleculeColumn, descriptors: descriptors});

    this.grid.dataFrame = this.template.enrichedDataFrame!;
  }
}