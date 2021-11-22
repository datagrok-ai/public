import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { NglAspect } from "./ngl-aspect"
import { PvizAspect } from "./pviz-aspect"
import { MiscMethods } from "./misc.js"

import { _package } from "../package";

export class TwinPviewer {
  root: HTMLElement;
  nglHost: HTMLElement;
  pVizHostL: HTMLElement;
  pVizHostH: HTMLElement;

  //representations: string[];
  repChoice: DG.InputBase;
  cdrScheme: DG.InputBase;
  ptmChoices: DG.InputBase;
  ptmMotifChoices: DG.InputBase;
  ptmProb: DG.InputBase;
  paratopes: DG.InputBase;
  sequenceTabs: DG.TabControl;

  panelNode: DG.DockNode;
  nglNode: DG.DockNode;
  sequenceNode: DG.DockNode;

  // ptmPredictions: string[];
  // ptmMotifPredictions: string[];
  twinSelections = { 'H': {}, 'L': {} };
  colorScheme = {
    "col_background": 'string',
    "col_heavy_chain": '#0069a7',
    "col_light_chain": '#f1532b',
    "col_cdr": '#45d145',
    "col_para": '#b0c4de',
    "col_highlight": '#45d145',
    "col_highlight_cdr": '#FFFF00',
    "col_partopes_low": '(176,196,222)', //col_para in rgb
    "col_partopes_high": '(255, 0, 255)'
  };

  openPanels;
  ngl;
  pViz;

  init(mlbView: DG.TableView, json, pdbStr: string) {
    // ---- SIDEPANEL REMOVAL ----
    let windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;

    // ---- INPUTS ----
    const representations = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    this.repChoice = ui.choiceInput('Representation', 'cartoon', representations);

    let schemes_lst = MiscMethods.extract_schemes(json);
    this.cdrScheme = ui.choiceInput('CDR3 Scheme', 'default', schemes_lst);

    let ptmKeys = [...new Set([...Object.keys(json.ptm_predictions.H), ...Object.keys(json.ptm_predictions.L)])];
    let ptmPredictions: string[] = [];
    let ptmMotifPredictions: string[] = [];

    for (let i = 0; i < ptmKeys.length; i++) {
      let ptmH = json.ptm_predictions.H[ptmKeys[i]];
      let ptmL = json.ptm_predictions.L[ptmKeys[i]];

      if ((typeof (ptmH) != "undefined" && ptmH[0][1] > 1) || (typeof (ptmL) != "undefined" && ptmL[0][1] > 1)) {
        ptmMotifPredictions.push(ptmKeys[i].replaceAll("_", " "));
      } else {
        ptmPredictions.push(ptmKeys[i].replaceAll("_", " "));
      }
    }

    //@ts-ignore
    this.ptmChoices = ui.multiChoiceInput('', [], ptmPredictions);
    //@ts-ignore
    this.ptmMotifChoices = ui.multiChoiceInput('', [], ptmMotifPredictions);
    this.ptmProb = ui.floatInput('PTM probability', 0.2);
    this.paratopes = ui.boolInput('Paratopes', false);

    // ---- INPUTS PANEL ----
    this.root = ui.div();
    let accOptions = ui.accordion();
    accOptions.addPane('3D model', () => ui.inputs([this.repChoice, this.cdrScheme]));
    accOptions.addPane('Sequence', () => ui.inputs([this.paratopes, this.ptmProb]));
    accOptions.addPane('Predicted PTMs', () => ui.div([this.ptmChoices]));
    accOptions.addPane('Motif PTMs', () => ui.div([this.ptmMotifChoices]));

    this.root.append(accOptions.root);

    // ---- VIEWER CONTAINERS ----
    this.nglHost = ui.div([], 'd4-ngl-viewer');
    this.pVizHostL = ui.box();
    this.pVizHostH = ui.box();

    //@ts-ignore
    this.sequenceTabs = ui.tabControl({
      'HEAVY': this.pVizHostH,
      'LIGHT': this.pVizHostL
    }).root;


    // ---- DOCKING ----
    this.panelNode = mlbView.dockManager.dock(this.root, 'right', null, 'NGL');
    this.nglNode = mlbView.dockManager.dock(this.nglHost, 'left', this.panelNode, 'NGL');
    //@ts-ignore
    this.sequenceNode = mlbView.dockManager.dock(this.sequenceTabs, 'down', this.nglNode, 'Sequence', 0.225);

    return [this.panelNode, this.nglNode, this.sequenceNode];

  }

  async reset(mlbView: DG.TableView) {
    if (!!this.ngl)
      this.ngl.stage.removeAllComponents();
    if (!!this.sequenceTabs)
    //@ts-ignore
      mlbView.dockManager.close(this.sequenceTabs);
    if (!!this.openPanels)
      this.openPanels.forEach((p) => mlbView.dockManager.close(p));
  }

  async twin(mlbView: DG.TableView, jsonStr, pdbStr) {

    this.openPanels = await this.init(mlbView, jsonStr, pdbStr);
    this.ngl = new NglAspect();
    await this.ngl.init(mlbView, this, pdbStr, jsonStr);
    let pViz = new PvizAspect();
    this.openPanels.push(await pViz.init(mlbView, this, this.ngl, jsonStr));

    this.repChoice.onChanged(async () => {
      await pViz.loadSequence(this, 'H', jsonStr, true)
      await pViz.loadSequence(this, 'L', jsonStr, true)
    });
    this.cdrScheme.onChanged(async () => {
      pViz.pVizParams.cdrMap = pViz.cdrMapping(this.cdrScheme.value, jsonStr)
      await pViz.loadSequence(this, 'H', jsonStr)
      await pViz.loadSequence(this, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    });
    this.paratopes.onChanged(async () => {
      await pViz.loadSequence(this, 'H', jsonStr)
      await pViz.loadSequence(this, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    });
    this.ptmChoices.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.ptmChoices.value, this.ptmProb.value, jsonStr)
      await pViz.loadSequence(this, 'H', jsonStr)
      await pViz.loadSequence(this, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    });
    this.ptmMotifChoices.onChanged(async () => {
      pViz.pVizParams.ptmMotifsMap = pViz.ptmMotifsMapping(this.ptmMotifChoices.value, this.ptmProb.value, jsonStr)
      await pViz.loadSequence(this, 'H', jsonStr)
      await pViz.loadSequence(this, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    });
    this.ptmProb.onChanged(async () => {
      pViz.pVizParams.ptmMap = pViz.ptmMapping(this.ptmChoices.value, this.ptmProb.value, jsonStr)
      await pViz.loadSequence(this, 'H', jsonStr)
      await pViz.loadSequence(this, 'L', jsonStr)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    });
  }
}