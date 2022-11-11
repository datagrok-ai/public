import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {NglAspect} from './ngl-aspect';
import {PvizAspect} from './pviz-aspect';
import {MiscMethods} from './misc';

export class TwinProteinView {
  root: HTMLElement;
  nglHost: HTMLElement;
  pVizHosts: { [key: string]: HTMLDivElement };
  repChoice: DG.InputBase;
  sequenceTabs: DG.TabControl;
  accOptions: DG.Accordion;
  panelNode: DG.DockNode;
  nglNode: DG.DockNode;
  sequenceNode: DG.DockNode;

  twinSelections: { [key: string]: {} };
  colorScheme = {
    'col_background': 'white',
    'col_chain': '#0069a7',
    'col_helix': 'red',
    'col_highlight': '#45d145',
  };

  ngl: NglAspect;
  pViz: PvizAspect;

  entry: string;
  ligandSelection: { [key: string]: any }

  public init(entry: string, bsView: DG.TableView, ligandSelection: { [key: string]: boolean }) {
    // ---- SIDEPANEL REMOVAL ----
    const windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;
    this.entry = entry;
    this.ligandSelection = ligandSelection;

    // ---- INPUTS ----
    const representations = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    this.repChoice = ui.choiceInput('Representation', 'cartoon', representations);
    this.root = ui.div();
    this.changeChoices();

    // ---- VIEWER CONTAINERS ----
    this.nglHost = ui.div([], 'd4-ngl-viewer');
    this.pVizHosts = {};
    this.twinSelections = {};


    // this.entry.entities[0].chains.forEach((chain) => {
    //   this.pVizHosts[chain.id] = ui.box();
    //   this.twinSelections[chain.id] = {};
    // });

    //@ts-ignore
    this.sequenceTabs = ui.tabControl(this.pVizHosts).root;

    this.ngl = new NglAspect();
    //this.pViz = new PvizAspect();

    this.panelNode = bsView.dockManager.dock(this.root, 'right', null, 'NGL');
    this.nglNode = bsView.dockManager.dock(this.nglHost, 'left', this.panelNode, 'NGL');
    //@ts-ignore
    this.sequenceNode = bsView.dockManager.dock(this.sequenceTabs, 'down', this.nglNode, 'Sequence', 0.225);
  }

  public async reset(entry: string) {
    this.entry = entry;
    this.twinSelections = {};
    this.pVizHosts = {};
    // entry.entities[0].chains.forEach((chain) => {
    //   this.twinSelections[chain.id] = {};
    //   this.pVizHosts[chain.id] = ui.box();
    // });

    const groups: { [_: string]: any } = {};
    const items: DG.TreeViewNode[] = [];

    for (const g in groups) {
      if (Object.prototype.hasOwnProperty.call(groups, g)) {
        groups[g].checked = false;
      }
    }
    for (const i of items) i.checked = false;

    this.changeChoices();

    if (!!this.ngl) {
      this.ngl.stage.removeAllComponents();
    }
  }

  public async show(bsView: DG.TableView) {
    const reload = async (val: boolean) => {
      // this.entry.entities[0].chains.forEach((chain) => {
      //   this.pViz.render(chain.id);
      // });

      this.ngl.render(val, this.ligandSelection);
      MiscMethods.setDockSize(bsView, this.nglNode, this.sequenceTabs);
    };

    await this.ngl.init(
      bsView,
      this.entry,
      this.colorScheme,
      this.nglHost,
      this.repChoice,
      this.twinSelections,
      this.ligandSelection,
    );
    //await this.pViz.init(this.entry, this.colorScheme, this.pVizHosts, this.twinSelections);

    MiscMethods.setDockSize(bsView, this.nglNode, this.sequenceTabs);

    grok.events.onCustomEvent('selectionChanged').subscribe((v) => {
      this.ngl.render(false, this.ligandSelection);
    });

    this.repChoice.onChanged(async () => {
      this.ngl.repChoice = this.repChoice;
      reload(true);
    });
  }

  public async changeLigands(bsView: DG.TableView, ligandSelection: { [key: string]: boolean }) {
    ligandSelection = ligandSelection;
    const reload = async (val: boolean) => {
      this.ngl.render(val, ligandSelection);
      MiscMethods.setDockSize(bsView, this.nglNode, this.sequenceTabs);
    };

    reload(false);
  }

  private changeChoices(): void {
    // ---- INPUTS PANEL ----
    if (!this.accOptions) {
      this.accOptions = ui.accordion();
    } else {
      this.accOptions.removePane(this.accOptions.getPane('3D model'));
    }

    this.accOptions.addPane('3D model', () => ui.inputs([this.repChoice]));
    this.root.append(this.accOptions.root);
  }
}
