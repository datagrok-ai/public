import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {NglAspect} from './ngl-aspect';
import {PvizAspect} from './pviz-aspect';
import {MiscMethods} from './misc';
import {ColorSchemeType, DataLoader, JsonType, ObsPtmType, PdbType} from '../utils/data-loader';
import {MlbEvents} from '../const';
import {Subscription} from 'rxjs';

export class TwinPviewer {
  dataLoader: DataLoader;

  root: HTMLElement;
  nglHost: HTMLElement;
  pVizHostL: HTMLDivElement;
  pVizHostH: HTMLDivElement;

  //representations: string[];
  repChoiceInput: DG.InputBase;
  cdrSchemeInput: DG.InputBase;
  ptmChoicesInput: DG.InputBase;
  ptmMotifChoicesInput: DG.InputBase;
  ptmObsChoicesInput: DG.InputBase;
  ptmProbInput: DG.InputBase;
  paratopesInput: DG.InputBase;
  sequenceTabs: DG.TabControl;

  schemesList: string[];

  accOptions: DG.Accordion;
  panelNode: DG.DockNode;
  nglNode: DG.DockNode;
  sequenceNode: DG.DockNode;

  // ptmPredictions: string[];
  // ptmMotifPredictions: string[];
  twinSelections = {'H': {}, 'L': {}};
  colorScheme: ColorSchemeType = {
    'colBackground': 'white',
    'colHeavyChain': '#0069a7',
    'colLightChain': '#f1532b',
    'colCdr': '#45d145',
    'colPara': '#b0c4de',
    'colHighlight': '#45d145',
    'colHighlightCdr': '#FFFF00',
    'colParatopesLow': 'rgb(176,196,222)', //col_para in rgb
    'colParatopesHigh': 'rgb(255,0,255)',
  };

  ngl: NglAspect;
  pViz: PvizAspect;

  isOpened: boolean = false;
  jsonStr: JsonType; // TODO: descriptive name
  pdbStr: PdbType;
  jsonObsPtm: ObsPtmType;

  subs: Subscription[] = [];

  constructor(dataLoader: DataLoader) {
    this.dataLoader = dataLoader;
  }

  public init(json: JsonType, pdb: PdbType, jsonObsPtm: ObsPtmType) {
    // ---- SIDEPANEL REMOVAL ----
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;

    this.jsonStr = json;
    this.pdbStr = pdb;
    this.jsonObsPtm = jsonObsPtm;

    // ---- INPUTS ----
    const representations = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    this.repChoiceInput = ui.choiceInput('Representation', 'cartoon', representations);

    this.schemesList = MiscMethods.extractSchemes(json);
    this.cdrSchemeInput = ui.choiceInput('CDR3 Scheme', MiscMethods.NoSchemeItem, this.schemesList);

    this.root = ui.div();
    this.changeChoices();

    this.ptmProbInput = ui.floatInput('PTM probability', 0.2);
    this.paratopesInput = ui.boolInput('Paratopes', false);

    // ---- VIEWER CONTAINERS ----
    this.nglHost = ui.div([], 'd4-ngl-viewer');
    this.pVizHostL = ui.box();
    this.pVizHostH = ui.box();

    this.sequenceTabs = ui.tabControl({
      'HEAVY': this.pVizHostH,
      'LIGHT': this.pVizHostL,
    });

    this.ngl = new NglAspect();
    this.pViz = new PvizAspect(this.dataLoader);
  }

  public async reset(json: JsonType, pdb: PdbType, jsonObsPtm: ObsPtmType) {
    this.jsonStr = json;
    this.pdbStr = pdb;
    this.jsonObsPtm = jsonObsPtm;

    this.twinSelections = {'H': {}, 'L': {}};

    const groups: { [_: string]: any } = {};
    const items: DG.TreeViewNode[] = [];

    for (const g in groups)
      if (groups[g].checked) groups[g].checked = false;

    for (const i of items) i.checked = false;

    this.changeChoices();

    if (!!this.ngl)
      this.ngl.stage.removeAllComponents();
  }

  public async open(mlbView: DG.TableView) {
    // ---- DOCKING ----
    if (!this.isOpened) {
      this.isOpened = true;
      this.panelNode = mlbView.dockManager.dock(this.root, DG.DOCK_TYPE.RIGHT, null, 'NGL', 0.2);
      this.nglNode = mlbView.dockManager.dock(this.nglHost, DG.DOCK_TYPE.LEFT, this.panelNode, 'NGL', 0.3);
      this.sequenceNode =
        mlbView.dockManager.dock(this.sequenceTabs.root, DG.DOCK_TYPE.DOWN, this.nglNode, 'Sequence', 0.225);
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs.root);

      this.subs.push(grok.events.onCustomEvent(MlbEvents.CdrChanged).subscribe((value: string) => {
        const key = this.schemesList.find((v) => v.toUpperCase() == value.toUpperCase());
        console.debug(`MLB: CompositionPviewer.onCustomEvent(${MlbEvents.CdrChanged}) ` +
          `value="${value}" -> key="${key}".`);
        this.cdrSchemeInput.value = key ? key : MiscMethods.NoSchemeItem;
      }));
    }
  }

  public async close(mlbView: DG.TableView) {
    if (!!this.sequenceTabs)
      mlbView.dockManager.close(this.sequenceTabs.root);

    if (!!this.panelNode)
      mlbView.dockManager.close(this.panelNode);

    if (!!this.nglNode)
      mlbView.dockManager.close(this.nglNode);

    if (!!this.sequenceNode)
      mlbView.dockManager.close(this.sequenceNode);

    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];

    this.isOpened = false;
  }

  public async show(mlbView: DG.TableView) {
    const doReload = async (reload: boolean) => {
      this.pViz.render('H');
      this.pViz.render('L');
      this.ngl.render(reload);
      MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs.root);
    };

    await this.ngl.init(mlbView, this.pdbStr, this.jsonStr, this.colorScheme, this.nglHost,
      this.repChoiceInput, this.cdrSchemeInput, this.paratopesInput, this.twinSelections);

    const obsChoice = this.ptmObsChoicesInput !== undefined ? this.ptmObsChoicesInput.value : null;

    await this.pViz.init(this.jsonStr, this.jsonObsPtm, this.colorScheme, this.pVizHostH, this.pVizHostL,
      this.ptmChoicesInput.value, this.ptmMotifChoicesInput.value, obsChoice, this.cdrSchemeInput,
      this.paratopesInput, this.ptmProbInput.value, this.twinSelections);

    MiscMethods.setDockSize(mlbView, this.nglNode, this.sequenceTabs.root);

    grok.events.onCustomEvent('selectionChanged').subscribe((v) => {
      this.ngl.render(false);
    });

    this.repChoiceInput.onChanged(async () => {
      this.ngl.repChoice = this.repChoiceInput;
      doReload(true);
    });

    this.cdrSchemeInput.onChanged(async () => {
      this.ngl.cdrScheme = this.cdrSchemeInput;
      this.pViz.cdrMapping(this.cdrSchemeInput);
      doReload(false);
    });

    this.paratopesInput.onChanged(async () => {
      this.ngl.paratopes = this.paratopesInput;
      this.pViz.parMapping(this.paratopesInput);
      doReload(false);
    });

    this.ptmChoicesInput.onChanged(async () => {
      this.pViz.ptmMapping(this.ptmChoicesInput.value, this.ptmProbInput.value);
      doReload(false);
    });

    this.ptmProbInput.onChanged(async () => {
      this.pViz.ptmMapping(this.ptmChoicesInput.value, this.ptmProbInput.value);
      doReload(false);
    });

    this.ptmMotifChoicesInput.onChanged(async () => {
      this.pViz.motMapping(this.ptmMotifChoicesInput.value, this.ptmProbInput.value);
      doReload(false);
    });

    if (this.ptmObsChoicesInput !== undefined) {
      this.ptmObsChoicesInput.onChanged(async () => {
        this.pViz.obsMapping(this.ptmObsChoicesInput.value);
        doReload(false);
      });
    }
  }

  private changeChoices(): void {
    const ptmKeys =
      [...new Set([...Object.keys(this.jsonStr.ptm_predictions.H), ...Object.keys(this.jsonStr.ptm_predictions.L)])];
    const ptmPredictions: string[] = [];
    const ptmMotifPredictions: string[] = [];

    for (let i = 0; i < ptmKeys.length; i++) {
      const ptmH = this.jsonStr.ptm_predictions.H[ptmKeys[i]];
      const ptmL = this.jsonStr.ptm_predictions.L[ptmKeys[i]];

      if ((typeof (ptmH) != 'undefined' && ptmH[0][1] > 1) || (typeof (ptmL) != 'undefined' && ptmL[0][1] > 1))
        ptmMotifPredictions.push(ptmKeys[i].replaceAll('_', ' '));
      else
        ptmPredictions.push(ptmKeys[i].replaceAll('_', ' '));
    }

    //@ts-ignore
    this.ptmChoicesInput = ui.multiChoiceInput('', [], ptmPredictions);
    //@ts-ignore
    this.ptmMotifChoicesInput = ui.multiChoiceInput('', [], ptmMotifPredictions);

    // ---- INPUTS PANEL ----
    if (!this.accOptions) {
      this.accOptions = ui.accordion();
    } else {
      this.accOptions.removePane(this.accOptions.getPane('3D model'));
      this.accOptions.removePane(this.accOptions.getPane('Sequence'));
      this.accOptions.removePane(this.accOptions.getPane('Predicted PTMs'));
      this.accOptions.removePane(this.accOptions.getPane('Motif PTMs'));

      if (this.accOptions.getPane('Observed PTMs') !== undefined)
        this.accOptions.removePane(this.accOptions.getPane('Observed PTMs'));
    }

    this.accOptions.addPane('3D model', () => ui.inputs([this.repChoiceInput, this.cdrSchemeInput]));
    this.accOptions.addPane('Sequence', () => ui.inputs([this.paratopesInput, this.ptmProbInput]));
    this.accOptions.addPane('Predicted PTMs', () => ui.div([this.ptmChoicesInput]));
    this.accOptions.addPane('Motif PTMs', () => ui.div([this.ptmMotifChoicesInput]));

    if (this.jsonObsPtm !== null) {
      const obsPtmKeys =
        [...new Set([...Object.keys(this.jsonObsPtm.H), ...Object.keys(this.jsonObsPtm.L)])];
      const obsPtmPredictions: string[] = [];

      for (let i = 0; i < obsPtmKeys.length; i++)
        obsPtmPredictions.push(obsPtmKeys[i].replaceAll('_', ' '));

      //@ts-ignore
      this.ptmObsChoicesInput = ui.multiChoiceInput('', [], obsPtmPredictions);
      this.accOptions.addPane('Observed PTMs', () => ui.div([this.ptmObsChoicesInput]));
    }

    this.root.append(this.accOptions.root);
  }
}
