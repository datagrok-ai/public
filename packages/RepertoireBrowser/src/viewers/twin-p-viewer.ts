import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { MolecularLiabilityBrowserPanel } from "./molecular-liability-browser-panel.js"
import { NglMethods } from "./ngl.js"
import { PvizMethods } from "./pViz.js"
import { MiscMethods } from "./misc.js"


import { _package } from "../package";

export class TwinPviewer {

  inputs;
  openPanels;
  ngl;
  pViz;

  async reset(mlbView: DG.TableView) {
    if (!!this.ngl)
      this.ngl.stage.removeAllComponents();
    if (!!this.inputs && !!this.inputs.sequence_tabs)
      mlbView.dockManager.close(this.inputs.sequence_tabs);
    if (!!this.openPanels)
      this.openPanels.forEach((p) => mlbView.dockManager.close(p));
  }

  async twin(mlbView: DG.TableView, jsonStr, pdbStr) {

    this.inputs = new MolecularLiabilityBrowserPanel();
    this.openPanels = await this.inputs.init(mlbView, jsonStr);
    this.ngl = new NglMethods();
    await this.ngl.init(mlbView, this.inputs, pdbStr, jsonStr);
    let pViz = new PvizMethods();
    this.openPanels.push(await pViz.init(mlbView, this.inputs, this.ngl, jsonStr));


    this.inputs.repChoice.onChanged(async () => {
      await pViz.loadSequence(this.inputs, 'H', jsonStr, true)
      await pViz.loadSequence(this.inputs, 'L', jsonStr, true)
    });
    this.inputs.cdr_scheme.onChanged(async () => {
      pViz.pVizParams.cdrMap = pViz.cdrMapping(this.inputs.cdr_scheme.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.paratopes.onChanged(async () => {
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_choices.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_motif_choices.onChanged(async () => {
      pViz.pVizParams.ptmMotifsMap = pViz.ptmMotifsMapping(this.inputs.ptm_motif_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
    this.inputs.ptm_prob.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.inputs.ptm_choices.value, this.inputs.ptm_prob.value, jsonStr)
      await pViz.loadSequence(this.inputs, 'H', jsonStr)
      await pViz.loadSequence(this.inputs, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.inputs.ngl_node, this.inputs.sequence_tabs);
    });
  }
}