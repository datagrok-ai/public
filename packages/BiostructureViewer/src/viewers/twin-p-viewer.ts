import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";

import { NglAspect } from "./ngl-aspect"
import { PvizAspect } from "./pviz-aspect"
import { MiscMethods } from "./misc"
import { PdbEntry } from '../pdb-entry';

import { _package } from "../package";


export class TwinPviewer {
  root: HTMLElement;
  nglHost: HTMLElement;
  pVizHosts: {[key: string]: HTMLDivElement};
  repChoice: DG.InputBase;
  sequenceTabs: DG.TabControl;
  accOptions: DG.Accordion;
  panelNode: DG.DockNode;
  nglNode: DG.DockNode;
  sequenceNode: DG.DockNode;

  twinSelections: {[key: string]: {}};
  colorScheme = {
    "col_background": 'white',
    "col_chain": '#0069a7',
    "col_helix": 'red',
    "col_highlight": '#45d145'
  };

  ngl: NglAspect;
  pViz: PvizAspect;

  entry: PdbEntry;

  public init(entry: PdbEntry, bsView: DG.TableView) {
    // ---- SIDEPANEL REMOVAL ----
    let windows = grok.shell.windows;
    windows.showProperties = false;
    windows.showHelp = false;
    windows.showConsole = false;
    this.entry = entry;

    // ---- INPUTS ----
    const representations = ['cartoon', 'backbone', 'ball+stick', 'licorice', 'hyperball', 'surface'];
    this.repChoice = ui.choiceInput('Representation', 'cartoon', representations);
    this.root = ui.div();
    this.changeChoices();

    // ---- VIEWER CONTAINERS ----
    this.nglHost = ui.div([], 'd4-ngl-viewer');
    this.pVizHosts = {};
    this.twinSelections = {};
    this.entry.entities[0].chains.forEach(chain => {
      this.pVizHosts[chain.id] = ui.box();
      this.twinSelections[chain.id] = {};
    });

    //@ts-ignore
    this.sequenceTabs = ui.tabControl(this.pVizHosts).root;

    this.ngl = new NglAspect();
    this.pViz = new PvizAspect();

    this.panelNode = bsView.dockManager.dock(this.root, 'right', null, 'NGL');
    this.nglNode = bsView.dockManager.dock(this.nglHost, 'left', this.panelNode, 'NGL');
    //@ts-ignore
    this.sequenceNode = bsView.dockManager.dock(this.sequenceTabs, 'down', this.nglNode, 'Sequence', 0.225);
  }

  public async reset(entry: PdbEntry) {
    this.entry = entry;
    this.twinSelections = {};
    this.pVizHosts = {};
    entry.entities[0].chains.forEach(chain => {
      this.twinSelections[chain.id] = {};
      this.pVizHosts[chain.id] = ui.box();
    });

    let groups: { [_: string]: any } = {};
    let items: DG.TreeViewNode[] = [];

    for (let g in groups) groups[g].checked = false;
    for (let i of items) i.checked = false;

    this.changeChoices();

    if (!!this.ngl)
      this.ngl.stage.removeAllComponents();
  }

  public async show(bsView: DG.TableView) {

    let reload = async (val: boolean) => {
      this.entry.entities[0].chains.forEach(chain => {
        this.pViz.render(chain.id);
      });

      this.ngl.render(val);
      MiscMethods.setDockSize(bsView, this.nglNode, this.sequenceTabs);
    };

    await this.ngl.init(bsView, this.entry, this.colorScheme, this.nglHost, this.repChoice, this.twinSelections);
    await this.pViz.init(this.entry, this.colorScheme, this.pVizHosts, this.twinSelections);

    MiscMethods.setDockSize(bsView, this.nglNode, this.sequenceTabs);

    grok.events.onCustomEvent("selectionChanged").subscribe((v) => { this.ngl.render(false); });

    this.repChoice.onChanged(async () => {
      this.ngl.repChoice = this.repChoice;
      reload(true);
    });

    // this.cdrScheme.onChanged(async () => {
    //   this.ngl.cdrScheme = this.cdrScheme;
    //   this.pViz.cdrMapping(this.cdrScheme);
    //   reload(false);
    // });
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