import {HitTriageBaseView} from "./hit-triage-base-view";
import * as DG from "datagrok-api/dg";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import {HitTriageApp} from "../hit-triage-app";


/**
 * Enrichment of the molecular dataset.
 **/
export class ComputeView extends HitTriageBaseView {

  grid: DG.Grid;

  constructor(app: HitTriageApp) {
    super(app);
    this.name = 'Hit Triage | Compute';

    this.grid = this.template.hitsTable!.plot.grid()
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
    this.template.enrichedTable = this.template.hitsTable;

    await this.template.enrichedTable!.columns.addNewCalculated("length", "length(${smiles})");

    // descriptors
    let descriptors = ["HeavyAtomCount", "NHOHCount"];
    await grok.chem.descriptors(this.template.enrichedTable!, this.template.hitsMolColumn, descriptors);

    // another way to invoke calculations
    //await grok.functions.call('ChemDescriptors', {table: session.enrichedDataFrame, smiles: session.sourceMoleculeColumn, descriptors: descriptors});

    this.grid.dataFrame = this.template.enrichedTable!;
  }
}