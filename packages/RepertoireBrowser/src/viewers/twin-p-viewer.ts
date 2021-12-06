import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { NglAspect } from "./ngl-aspect"
import { PvizAspect } from "./pviz-aspect"
import { MiscMethods } from "./misc.js"

import { _package } from "../package";
import { resolveModuleName } from "typescript";

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
  ptmObsChoices: DG.InputBase;
  ptmProb: DG.InputBase;
  paratopes: DG.InputBase;
  sequenceTabs: DG.TabControl;

  accOptions: DG.Accordion;
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

  ngl: NglAspect;
  pViz: PvizAspect;

  isOpened: boolean = false;
  jsonStr: any;
  pdbStr: string;
  jsonObsPtm: any;

  init(json: any, pdb: string, jsonObsPtm: any) {
    // ---- SIDEPANEL REMOVAL ----
    let windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;

    this.jsonStr = json;
    this.pdbStr = pdb;
    this.jsonObsPtm = jsonObsPtm;

    // ---- INPUTS ----
    const representations = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    this.repChoice = ui.choiceInput('Representation', 'cartoon', representations);

    let schemes_lst = MiscMethods.extract_schemes(json);
    this.cdrScheme = ui.choiceInput('CDR3 Scheme', 'default', schemes_lst);

    this.root = ui.div();
    this.changeChoices();

    this.ptmProb = ui.floatInput('PTM probability', 0.2);
    this.paratopes = ui.boolInput('Paratopes', false);

    // ---- VIEWER CONTAINERS ----
    this.nglHost = ui.div([], 'd4-ngl-viewer');
    this.pVizHostL = ui.box();
    this.pVizHostH = ui.box();

    //@ts-ignore
    this.sequenceTabs = ui.tabControl({
      'HEAVY': this.pVizHostH,
      'LIGHT': this.pVizHostL
    }).root;

    this.ngl = new NglAspect();
    this.pViz = new PvizAspect();
  }

  async reset(json: string, pdb: string, jsonObsPtm: any) {
    this.jsonStr = json;
    this.pdbStr = pdb;
    this.jsonObsPtm = jsonObsPtm;

    this.twinSelections = { 'H': {}, 'L': {} };

    let groups: { [_: string]: any } = {};
    let items: DG.TreeViewNode[] = [];

    for (let g in groups) groups[g].checked = false;
    for (let i of items) i.checked = false;

    this.changeChoices();

    if (!!this.ngl)
      this.ngl.stage.removeAllComponents();
  }

  async open(mlbView: DG.TableView) {
    // ---- DOCKING ----
    if (!this.isOpened) {
      this.isOpened = true;
      this.panelNode = mlbView.dockManager.dock(this.root, 'right', null, 'NGL');
      this.nglNode = mlbView.dockManager.dock(this.nglHost, 'left', this.panelNode, 'NGL');
      //@ts-ignore
      this.sequenceNode = mlbView.dockManager.dock(this.sequenceTabs, 'down', this.nglNode, 'Sequence', 0.225);
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    }
  }

  async close(mlbView: DG.TableView) {

    if (!!this.sequenceTabs)
      //@ts-ignore
      mlbView.dockManager.close(this.sequenceTabs);

    if (!!this.panelNode)
      mlbView.dockManager.close(this.panelNode);

    if (!!this.nglNode)
      mlbView.dockManager.close(this.nglNode);

    if (!!this.sequenceNode)
      mlbView.dockManager.close(this.sequenceNode);

    this.isOpened = false;
  }

  async twin(mlbView: DG.TableView) {

    let reload = async(val: boolean) =>{
      await this.pViz.loadSequence(this, 'H', this.jsonStr, this.jsonObsPtm, val);
      await this.pViz.loadSequence(this, 'L', this.jsonStr, this.jsonObsPtm, val)
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs);
    };

    await this.ngl.init(mlbView, this, this.pdbStr, this.jsonStr);
    await this.pViz.init(mlbView, this, this.ngl, this.jsonStr, this.jsonObsPtm);

    this.repChoice.onChanged(async () => {
      reload(true);
    });
    this.cdrScheme.onChanged(async () => {
      this.pViz.pVizParams.cdrMap = this.pViz.cdrMapping(this.cdrScheme.value, this.jsonStr)
      reload(false);
    });
    this.paratopes.onChanged(async () => {
      reload(false);
    });
    this.ptmChoices.onChanged(async () => {
      this.pViz.pVizParams.ptmMap = this.pViz.ptmMapping(this.ptmChoices.value, this.ptmProb.value, this.jsonStr)
      reload(false);
    });
    this.ptmMotifChoices.onChanged(async () => {
      this.pViz.pVizParams.ptmMotifsMap = this.pViz.ptmMotifsMapping(this.ptmMotifChoices.value, this.ptmProb.value, this.jsonStr)
      reload(false);
    });

    if(this.ptmObsChoices !== undefined){
      this.ptmObsChoices.onChanged(async () => {
        this.pViz.pVizParams.ptmObsMap = this.pViz.ptmObsMapping(this.ptmObsChoices.value, this.jsonObsPtm)
        reload(false);
      });
    }

    this.ptmProb.onChanged(async () => {
      this.pViz.pVizParams.ptmMap = this.pViz.ptmMapping(this.ptmChoices.value, this.ptmProb.value, this.jsonStr)
      reload(false);
    });
  }

  changeChoices(): void {

    let ptmKeys = [...new Set([...Object.keys(this.jsonStr.ptm_predictions.H), ...Object.keys(this.jsonStr.ptm_predictions.L)])];
    let ptmPredictions: string[] = [];
    let ptmMotifPredictions: string[] = [];

    for (let i = 0; i < ptmKeys.length; i++) {
      let ptmH = this.jsonStr.ptm_predictions.H[ptmKeys[i]];
      let ptmL = this.jsonStr.ptm_predictions.L[ptmKeys[i]];

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

    // ---- INPUTS PANEL ----
    if (!this.accOptions) {
      this.accOptions = ui.accordion();
    } else {
      this.accOptions.removePane(this.accOptions.getPane('3D model'));
      this.accOptions.removePane(this.accOptions.getPane('Sequence'));
      this.accOptions.removePane(this.accOptions.getPane('Predicted PTMs'));
      this.accOptions.removePane(this.accOptions.getPane('Motif PTMs'));

      if(this.accOptions.getPane('Observed PTMs') !== undefined)
        this.accOptions.removePane(this.accOptions.getPane('Observed PTMs'));
    }

    this.accOptions.addPane('3D model', () => ui.inputs([this.repChoice, this.cdrScheme]));
    this.accOptions.addPane('Sequence', () => ui.inputs([this.paratopes, this.ptmProb]));
    this.accOptions.addPane('Predicted PTMs', () => ui.div([this.ptmChoices]));
    this.accOptions.addPane('Motif PTMs', () => ui.div([this.ptmMotifChoices]));

    if (this.jsonObsPtm !== null) {
      let obsPtmKeys = [...new Set([...Object.keys(this.jsonObsPtm.ptm_observed.H), ...Object.keys(this.jsonObsPtm.ptm_observed.L)])];
      let obsPtmPredictions: string[] = [];

      for (let i = 0; i < obsPtmKeys.length; i++) {
        obsPtmPredictions.push(obsPtmKeys[i].replaceAll("_", " "));
      }
      //@ts-ignore
      this.ptmObsChoices = ui.multiChoiceInput('', [], obsPtmPredictions);
      this.accOptions.addPane('Observed PTMs', () => ui.div([this.ptmObsChoices]));
    } else {

    }

    this.root.append(this.accOptions.root);
  }
}