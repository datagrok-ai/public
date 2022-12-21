import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import $ from "cash-dom";
import {_rdKitModule} from '../utils/chem-common-rdkit';
import {isMolBlock} from "../utils/convert-notation-utils";
import {chem} from "datagrok-api/grok";
import {TreeViewGroup, TreeViewNode} from "datagrok-api/dg";
import Sketcher = chem.Sketcher;
import {RDKitCellRenderer} from "../rendering/rdkit-cell-renderer";
import {chemSubstructureSearchLibrary} from "../chem-searches";
import {getScaffoldTree, installScaffoldGraph} from "../package";

let renderer : RDKitCellRenderer | null= null;

function debugParent(node: TreeViewNode) : void {
  let count = 0;
  let parent = node.parent;
  for (; parent !== null && parent !== undefined && count < 100; parent = parent.parent, ++count) {
    console.log('DEBUG NODE PARENTS: ' + count + " " + parent.constructor.name + " " + ((parent.value as any).value !== undefined ? (parent.value as any).value.smiles : "NA"));
  }
}


function enableNodeExtendArrow(group: TreeViewGroup, enable: boolean): void {
  const c = group.root.getElementsByClassName('d4-tree-view-tri');
  if (c.length > 0)
    (c[0] as HTMLElement).style.visibility = enable ? 'visible' : 'hidden';
}

function buildOrphans(rootGroup: TreeViewGroup) {
  if (rootGroup.parent !== null && rootGroup.children.length > 0) { //folder groups don't have children, will skip

    const v = rootGroup.value as any;
    v.bitset_orphans = null;
    const bitsetThis = v.bitset;

    let bitsetParent = null;
    if(!rootGroup.parent || !rootGroup.parent.value) {
      bitsetParent = DG.BitSet.create(bitsetThis.length);
      bitsetParent = bitsetParent.or(bitsetThis, false);
    } else bitsetParent =  (rootGroup.parent.value as any).bitset;

    let bitsetOrphans = DG.BitSet.create(bitsetThis.length);
    bitsetOrphans = bitsetOrphans.or(bitsetThis, false);
    bitsetOrphans = bitsetOrphans.and(bitsetParent); //all matched in the parent

    let bitsetChild = null;
    for (let n = 0; n < rootGroup.children.length; ++n) {
      bitsetChild = (rootGroup.children[n].value as any).bitset;
      bitsetOrphans = bitsetOrphans.xor(bitsetChild, false);

      let bitsetMask = DG.BitSet.create(bitsetOrphans.length);
      bitsetMask = bitsetMask.or(bitsetOrphans, false);
      bitsetMask = bitsetMask.and(bitsetChild, false);
      bitsetOrphans = bitsetOrphans.xor(bitsetMask, false);
    }
    if (bitsetOrphans.trueCount > 0)
      v.bitset_orphans = bitsetOrphans;
  }

  for (let n = 0; n < rootGroup.children.length; ++n) {
    buildOrphans(rootGroup.children[n] as TreeViewGroup);
  }
}

function fillVisibleNodes(rootGroup: TreeViewGroup, visibleNodes: Array<TreeViewGroup>) : void {
  if (rootGroup.parent !== null) {
    visibleNodes.push(rootGroup);
    if (!rootGroup.expanded) return;
  }
  for (let n = 0; n < rootGroup.children.length; ++n) {
    fillVisibleNodes(rootGroup.children[n] as TreeViewGroup, visibleNodes);
  }
}

function updateNodeHitsLabel(group : TreeViewNode, text : string) : void {
  const labelDiv = (group.value as any).labelDiv;
  labelDiv!.innerHTML = text;
}

async function updateAllNodesHits(thisViewer: ScaffoldTreeViewer) {
  const items = thisViewer.tree.items;
  if (items.length > 0) {
    await updateNodesHitsImpl(thisViewer, items, 0, items.length - 1);
    console.log('All Hits are calculated');
    buildOrphans( thisViewer.tree);
    console.log('Orphans are calculated');
    thisViewer.appendOrphanFolders(thisViewer.tree);
    thisViewer.updateSizes();
    console.log('Orphans are appended');
  }
}

async function updateVisibleNodesHits(thisViewer: ScaffoldTreeViewer) {
  const visibleNodes : Array<TreeViewGroup> = [];
  fillVisibleNodes(thisViewer.tree, visibleNodes);

  const start = Math.floor(thisViewer.tree.root.scrollTop / CELL_HEIGHT);
  let end = start + Math.ceil(thisViewer.root.offsetHeight / CELL_HEIGHT);
  if (end >= visibleNodes.length)
    end = visibleNodes.length - 1;

  updateNodesHitsImpl(thisViewer, visibleNodes, start, end);
}

async function updateNodesHitsImpl(thisViewer: ScaffoldTreeViewer, visibleNodes : Array<TreeViewNode>, start: number, end: number) {
  let group : TreeViewNode | null = null;
  //console.log('scroll ' + start + ' ' + +end + ' total ' + visibleNodes.length);
  let v = null;
  for (let n = start; n <= end; ++n) {
    group = visibleNodes[n];
    v = group.value as any;
    if (v.init || v.orphans)
      continue;

    // @ts-ignore
    /* in case of delayed rendering
    const canvas: HTMLCanvasElement = (group.value as any).canvas;
    const r = window.devicePixelRatio;
    // @ts-ignore
    renderer.render(canvas.getContext('2d'), 0, 0, Math.floor(canvas.width / r), Math.floor(canvas.height / r), DG.GridCell.fromValue((group.value as any).smiles));
    */
    const bitset = thisViewer.molColumn === null ? null : await chemSubstructureSearchLibrary(thisViewer.molColumn,  v.smiles, '');
    v.bitset = bitset;
    v.init = true;
    updateNodeHitsLabel(group, bitset!.trueCount.toString());
  }
}

async function _initWorkers(molColumn: DG.Column) : Promise<DG.BitSet> {
  const molStr = molColumn.length === 0 ? '' : molColumn.get(molColumn.length -1);
  return await chemSubstructureSearchLibrary(molColumn, molStr, '');
}

function renderMolecule(molStr: string, width: number, height: number, skipDraw: boolean = false): HTMLDivElement {
  const moleculeHost = ui.canvas(width, height);
  $(moleculeHost).addClass('chem-canvas');
  const r = window.devicePixelRatio;
  moleculeHost.width = width;
  moleculeHost.height = height;
  moleculeHost.style.width = (width / r).toString() + 'px';
  moleculeHost.style.height= (height / r).toString() + 'px';
  if (renderer === null)
    renderer = new RDKitCellRenderer(_rdKitModule);

  if (skipDraw) {
    const g = moleculeHost.getContext('2d');
    g!.font = "18px Roboto, Roboto Local";
    const text = 'Loading...';
    let tm = g!.measureText(text);
    const fontHeight = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent
    const lineWidth = tm.width;
    g!.fillText(text, Math.floor((width - lineWidth) / 2), Math.floor((height - fontHeight) / 2));
  } else // @ts-ignore
   renderer.render(moleculeHost.getContext('2d')!, 0, 0, Math.floor(width / r), Math.floor(height / r), DG.GridCell.fromValue(molStr));

  return ui.divH([moleculeHost], 'chem-mol-box');
}


function getFlagIcon(group: TreeViewGroup) : HTMLElement | null {
  const molHost: HTMLElement = group.captionLabel;
  const c = molHost.getElementsByClassName('icon-fill');
 return c.length === 0 ? null : c[0] as HTMLElement;
}

const CELL_HEIGHT = 120;
const CELL_WIDTH = 200;
const GENERATE_ERROR_MSG = 'Generating tree failed...Please check the dataset';
const NO_MOL_COL_ERROR_MSG = 'There is no molecule column available';

export class ScaffoldTreeViewer extends DG.JsViewer {
  tree: DG.TreeViewGroup;
  bitset: DG.BitSet | null = null;
  wrapper: SketcherDialogWrapper | null = null;
  molColumns: Array<DG.Column> = [];
  molColumnIdx: number = -1;
  checkBoxesUpdateInProgress: boolean = false;
  TreeEncodeUpdateInProgress: boolean = false;
  _generateLink?: HTMLElement;
  _message?: HTMLElement | null = null;
  _iconAdd: HTMLElement | null = null;
  _iconDelete: HTMLElement | null = null;

  dirtyFlag: boolean = false;

  workersInit: boolean = false;

  MoleculeColumn: string;
  TreeEncode: string;

  constructor() {
    super();

    this.tree = ui.tree();

    const dframe = grok.shell.tv.dataFrame;
    this.molColumns = dframe.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
    const molColNames = new Array(this.molColumns.length);
    for (let n = 0; n < molColNames.length; ++n) {
      molColNames[n] = this.molColumns[n].name;
    }
    this.molColumnIdx = this.molColumns.length > 0 ? 0 : -1;

    this.MoleculeColumn = this.string('MoleculeColumn', molColNames.length === 0 ? null : molColNames[0], {
      choices: molColNames,
      category: 'Data',
      userEditable: this.molColumns.length > 0
    });
    this.TreeEncode = this.string('TreeEncode', '[]', {category: 'Data', userEditable: false});
  }

  get treeRoot() {
    return this.tree
  };

  private get message(): string {
    return this._message?.innerHTML as string;
  }

  private set message(msg: string | null) {
    if (this._message === undefined || this._message === null)
      return;

    this._message.style.visibility = msg === null ? 'hidden' : 'visible';
    // @ts-ignore
    this._message.innerHTML = msg ?? '';
  }

  async generateTree() {
    if (this.molColumn === null)
      return;

    this.dirtyFlag = false;
    this.message = null;
    this.clear();

    ui.setUpdateIndicator(this.root, true);
    const progressBar = DG.TaskBarProgressIndicator.create('Generating Scaffold Tree...');

    try {
      await installScaffoldGraph();
    } catch (e) {
      console.error(e);
      ui.setUpdateIndicator(this.root, false);
      progressBar.update(100, 'Installation of ScaffoldGraph failed');
      this.message = GENERATE_ERROR_MSG;
      //new Balloon().error('Could not generate scaffold tree');
      return;
    }

    if (!this.workersInit) {
      await _initWorkers(this.molColumn);
      this.workersInit = true;
    }

    const maxMolCount = 750;
    let length = this.molColumn.length;
    let step = 1;
    if (this.molColumn.length > maxMolCount) {
      step = Math.floor(this.molColumn.length / maxMolCount);
      if (step === 0) step = 1;
      length = maxMolCount;
    }

    const ar = new Array(length);
    for (let n = 0, m = 0; n < length; ++n, m += step) {
      ar[n] = this.molColumn.get(m);
    }

    const molCol: DG.Column = DG.Column.fromStrings('smiles', ar);
    molCol.semType = DG.SEMTYPE.MOLECULE;
    const dframe = DG.DataFrame.fromColumns([molCol]);

    let jsonStr = null;
    try {
      jsonStr = await getScaffoldTree(dframe);
    } catch (e) {
      console.error(e);
      ui.setUpdateIndicator(this.root, false);
      progressBar.update(50, 'Building Scaffold Tree Failed');
      this.message = GENERATE_ERROR_MSG;
      //new Balloon().error('Could not generate scaffold tree');
      return;
    }
    //progressBar.update(50, 'Initializing Scaffold Tree..: 80% completed');

    if (jsonStr !== null) {
      const json = JSON.parse(jsonStr);
      console.log(json);

      const thisViewer = this;
      ScaffoldTreeViewer.deserializeTrees(json, this.tree, (molStr: string, rootGroup: TreeViewGroup, countNodes: number) => {
        return thisViewer.createGroup(molStr, rootGroup, false, countNodes);
      });

      await updateVisibleNodesHits(this); //first visible N nodes
      ui.setUpdateIndicator(this.root, false);
      progressBar.update(100, 'Finished Building Scaffold Tree');

      this.updateSizes();
      this.updateUI();
      updateAllNodesHits(this); //this will run asynchronously
    }
  }

  get molColumn(): DG.Column | null {
    return this.molColumns.length === 0 ? null : this.molColumns[this.molColumnIdx];
  }

  private openEditSketcher(group: TreeViewGroup) {
    if (this.wrapper !== null) {
      this.wrapper.node = group;
      return;
    }

    //debugParent(group);

    const thisViewer = this;
    this.wrapper = SketcherDialogWrapper.create("Edit Scaffold...", "Save", group, async(molStrSketcher: string, node: TreeViewGroup, valid: boolean) => {

      while (node.captionLabel.firstChild) {
        node.captionLabel.removeChild(node.captionLabel.firstChild);
      }
      const bitset = thisViewer.molColumn === null ? null : await chemSubstructureSearchLibrary(thisViewer.molColumn, molStrSketcher, '');
      const molHost = renderMolecule(molStrSketcher, CELL_WIDTH, CELL_HEIGHT);
      this.addIcons(molHost, bitset!.trueCount.toString(), group);
      node.captionLabel.appendChild(molHost);
      const iconRoot = getFlagIcon(node);
      iconRoot!.style.cssText = iconRoot!.style.cssText += ('color: ' + (valid ? 'lightgreen !important' : 'hotpink !important'));
      iconRoot!.style.visibility = 'visible';

      node.value = {smiles: molStrSketcher, bitset: bitset};

      //orphans
      const parent = !node.parent ? node : node.parent as TreeViewGroup;
      buildOrphans(parent);
      thisViewer.appendOrphanFolders(parent);

      thisViewer.wrapper?.close();
      thisViewer.wrapper = null;
      thisViewer.updateFilters();
      thisViewer.dirtyFlag = true;
      thisViewer.TreeEncodeUpdateInProgress = true;
      thisViewer.TreeEncode = JSON.stringify(ScaffoldTreeViewer.serializeTrees(thisViewer.tree));
      thisViewer.TreeEncodeUpdateInProgress = false;
    }, (smilesSketcher: string, node: TreeViewGroup) => {

      if (node.parent === null)
        return true;

      let success = true;
      if (node.parent.value !== null) {
        const smilesParent = (node.parent.value as any).smiles;
        success = ScaffoldTreeViewer.validateNodes(smilesSketcher, smilesParent);
        if (!success)
          return false;
      }
      const children = node.items;
      for (let n = 0; n < children.length; ++n) {
        success = ScaffoldTreeViewer.validateNodes((children[n].value as any).smiles, smilesSketcher);
        if (!success)
          return false;
      }
      return true;
    }, () => {
      if (thisViewer.wrapper != null) {
        thisViewer.clearFilters();
        thisViewer.wrapper?.close();
        thisViewer.wrapper = null;
      }
    }, () => {
      if (thisViewer.wrapper != null) {
        thisViewer.clearFilters();
        thisViewer.wrapper?.close();
        thisViewer.wrapper = null;
      }
    }, async(strMolSketch: string) => {
      await thisViewer.filterByStruct(strMolSketch);
    });

    this.wrapper.show();
  }

  private openAddSketcher(group: TreeViewGroup) {
    const v = (group.value as any);
    const molStr = v === null ? "" : v.smiles;
    if (this.wrapper !== null) {
      this.wrapper.node = group;
      return;
    }
    const thisViewer = this;
    this.wrapper = SketcherDialogWrapper.create("Add New Scaffold...", "Add", group,async  (molStrSketcher: string, parent: TreeViewGroup, valid: boolean) => {
      debugParent(parent);
      const child = thisViewer.createGroup(molStrSketcher, parent);
      if (child !== null) {
        enableNodeExtendArrow(child, false);
        enableNodeExtendArrow(parent, true);
        const v = child.value as any;
        const bitset = thisViewer.molColumn === null ? null : await chemSubstructureSearchLibrary(thisViewer.molColumn, v.smiles, '');
        v.bitset = bitset;
        v.init = true;
        updateNodeHitsLabel(child, bitset!.trueCount.toString());
      }

      //orphans
      buildOrphans(parent);
      thisViewer.appendOrphanFolders(parent);

      thisViewer.updateSizes();
      thisViewer.updateUI();
      thisViewer.updateFilters();
      thisViewer.wrapper?.close();
      thisViewer.wrapper = null;
      thisViewer.dirtyFlag = true;
      thisViewer.TreeEncodeUpdateInProgress = true;
      thisViewer.TreeEncode = JSON.stringify(ScaffoldTreeViewer.serializeTrees(thisViewer.tree));
      thisViewer.TreeEncodeUpdateInProgress = false;
    }, (smilesSketcher: string) => {
      return ScaffoldTreeViewer.validateNodes(smilesSketcher, molStr);
    }, () => {
      if (thisViewer.wrapper != null) {
        thisViewer.clearFilters();
        thisViewer.wrapper?.close();
        thisViewer.wrapper = null;
      }
    }, () => {
      if (thisViewer.wrapper != null) {
        thisViewer.clearFilters();
        thisViewer.wrapper?.close();
        thisViewer.wrapper = null;
      }
    }, async(strMolSketch: string) => {
      await this.filterByStruct(strMolSketch);
    });
    this.wrapper.show();
  }

  clear() {
    this.clearFilters();
    while (this.tree.children.length > 0) {
      this.tree.children[0].remove();
    }

    this.TreeEncodeUpdateInProgress = true;
    this.TreeEncode = JSON.stringify(ScaffoldTreeViewer.serializeTrees(this.tree));
    this.TreeEncodeUpdateInProgress = false;

    this.updateUI();
  }

  clearFilters(): void {
    if (this.bitset === null)
      return;

    this.bitset = null;
    if (this.molColumn !== null)
      delete this.molColumn.temp['chem-scaffold-filter'];

    this.checkBoxesUpdateInProgress = true;
    const checkedNodes = this.tree.items.filter((v) => v.checked);
    for (let n = 0; n < checkedNodes.length; ++n)
      checkedNodes[n].checked = false;
    this.checkBoxesUpdateInProgress = false;

    this.dataFrame.rows.requestFilter();
    this.updateUI();
  }

  selectTableRows(group: TreeViewGroup, flag: boolean): void {
    const bitset = (group.value as any).bitset as DG.BitSet;
    if (flag)
      this.molColumn?.dataFrame.selection.or(bitset);
    else
      this.molColumn?.dataFrame.selection.andNot(bitset);
  }

  updateFilters(): void {
    if (this.molColumn === null)
      return;

    const checkedNodes = this.tree.items.filter((v) => v.checked);
    if(checkedNodes.length === 0) {
      this.clearFilters();
      return;
    }

    if (checkedNodes.length === 1) {
      const molStr = (checkedNodes[0].value as any).smiles;
      if(molStr !== undefined) {
        const mol = _rdKitModule.get_mol(molStr);
        const molFile = mol.get_molblock();
        mol.delete();
        this.molColumn.temp['chem-scaffold-filter'] = molFile;
      }
    } else delete this.molColumn.temp['chem-scaffold-filter'];

    if (this.bitset === null)
      this.bitset = DG.BitSet.create(this.molColumn.length);

    this.bitset.setAll(false, false);
    for (let n = 0; n < checkedNodes.length; ++n) {
      let bitset = (checkedNodes[n].value as any).bitset;
      this.bitset = this.bitset.or(bitset);
    }

    this.dataFrame.rows.requestFilter();
    this.updateUI();
  }

  private async filterByStruct(strMol: string) {
    if (this.molColumn != null && strMol !== null) {
      const mol = _rdKitModule.get_mol(strMol);
      const molFile = mol.get_molblock();
      mol.delete();
      this.molColumn.temp['chem-scaffold-filter'] = molFile;
      const bitset = await chemSubstructureSearchLibrary(this.molColumn, strMol, '');
      if (this.bitset === null)
        this.bitset = bitset;
      else {
        this.bitset.setAll(false, false);
        this.bitset = this.bitset.or(bitset);
      }

      this.dataFrame.rows.requestFilter();
      this.updateUI();
    }
  }

  addIcons(molHost: HTMLDivElement, label: string, group: TreeViewGroup): void {
    const iconsDiv = ui.divV([
      ui.iconFA('plus', () => this.openAddSketcher(group), 'Add new scaffold'),
      ui.iconFA('edit', () => this.openEditSketcher(group), 'Edit scaffold'),
      ui.divText(''),
      ui.iconFA('check-square', () => this.selectTableRows(group, true), 'Select rows'),
      ui.iconFA('square', () => this.selectTableRows(group, false), 'Unselect rows')
    ]);
    iconsDiv.onclick = (e) => e.stopImmediatePropagation();
    iconsDiv.onmousedown = (e) => e.stopImmediatePropagation();
    iconsDiv.style.justifyContent = 'center';

    const flagIcon = ui.iconFA('circle', () => this.openEditSketcher(group), 'The scaffold was edited');
    flagIcon.style.fontSize = '8px';
    flagIcon.style.marginLeft = '5px';
    //flagIcon.style.color = 'hotpink !important';
    flagIcon.style.visibility = 'hidden';
    flagIcon.classList.remove('fal');
    flagIcon.classList.add('fas', 'icon-fill');

    let labelDiv = null;
    const iconsInfo = ui.divH([labelDiv = ui.divText(label), flagIcon]);
    (group.value as any).labelDiv = labelDiv;

    const c = molHost.getElementsByTagName('CANVAS')
    if (c.length > 0)
      (group.value as any).canvas = c[0];

    iconsInfo.onclick = (e) => e.stopImmediatePropagation();
    iconsInfo.onmousedown = (e) => e.stopImmediatePropagation();
    molHost.appendChild(ui.divV([iconsInfo, iconsDiv]));
    ui.tools.setHoverVisibility(molHost, [iconsDiv]);
  }


  private createGroup(molStr: string, rootGroup: TreeViewGroup, skipDraw: boolean = false, index: number = -1) {
    if (this.molColumn === null)
      return null;

    const mol = _rdKitModule.get_mol(molStr);
    molStr = mol.get_molblock();
    mol.delete();

    const bitset =  DG.BitSet.create(this.molColumn.length);//await chemSubstructureSearchLibrary(this.molColumn, molStr, '');
    //const ratio = bitset.trueCount / bitset.length;
    //if (ratio !== 0 &&  bitset.trueCount === 1)//ratio < 0.05)
     // return null;
    const molHost = renderMolecule(molStr, CELL_WIDTH, CELL_HEIGHT, skipDraw);
    const group = rootGroup.group(molHost, {smiles: molStr, bitset: bitset, bitset_orphans : null});
    this.addIcons(molHost, bitset.trueCount.toString(), group);

    // Highlighting
    // molHost.onmouseenter = () => this.dataFrame.rows.match(bitset).highlight();
    // molHost.onmouseleave = () => this.dataFrame.rows.match('false').highlight();
    group.enableCheckBox(false);
    group.autoCheckChildren = false;
    const thisViewer = this;
    group.onNodeCheckBoxToggled.subscribe((node: TreeViewNode) => {
      if (!thisViewer.checkBoxesUpdateInProgress && node.value !== null)
        thisViewer.updateFilters();
    });

    return group;
  }

  createOrphansGroup(rootGroup: TreeViewGroup) {
    const divHost = ui.iconFA('folder');
    divHost.style.fontSize = '66px';
    divHost.style.width = '200px';
    divHost.style.height = '120px';
    divHost.style.top = '33px'; //0.5 from the height
    divHost.style.cssText += 'color: #50a9c573 !important';
    divHost.classList.remove('fal');
    divHost.classList.add('fas', 'icon-fill');

    const group = rootGroup.group(divHost, {orphans: true});
    if (group.children.length === 0)
      enableNodeExtendArrow(group, false);

    group.enableCheckBox(false);
    group.autoCheckChildren = false;
    const thisViewer = this;
    group.onNodeCheckBoxToggled.subscribe((node: TreeViewNode) => {
      if (!thisViewer.checkBoxesUpdateInProgress && node.value !== null)
        thisViewer.updateFilters();
    });

    return group;
  }

   appendOrphanFolders(rootGroup: TreeViewGroup) {
   // console.log(rootGroup.constructor.name);

    if (rootGroup.parent) {
      for (let n = 0; n < rootGroup.children.length; ++n) {
        if ((rootGroup.children[n].value as any).orphans)
          rootGroup.children[n].remove();
      }

      const bitsetOrphans = (rootGroup.value as any).bitset_orphans;
      if (bitsetOrphans !== null && bitsetOrphans !== undefined && bitsetOrphans.trueCount > 0) {
        const group = this.createOrphansGroup(rootGroup);
        (group.value as any).bitset = bitsetOrphans;
      }
    }
    for (let n = 0; n < rootGroup.children.length; ++n) {
      this.appendOrphanFolders(rootGroup.children[n] as TreeViewGroup);
    }
  }


  onPropertyChanged(p: DG.Property): void {
    if (p.name === 'MoleculeColumn') {
      this.clear();

      for (let n = 0; n < this.molColumns.length; ++n) {
        if (this.molColumns[n].name === this.MoleculeColumn) {
          this.molColumnIdx = n;
          break;
        }
      }

      //this.molColumn = this.dataFrame.columns.byName(this.MoleculeColumn);
      //console.log('Property changed: ' + p.name);
    } else if (p.name === 'TreeEncode') {
      if (this.TreeEncodeUpdateInProgress)
        return;

      const json = JSON.parse(this.TreeEncode);
      ScaffoldTreeViewer.deserializeTrees(json, this.tree, () => {
        (molStr: string, rootGroup: TreeViewGroup) => {
          return this.createGroup(molStr, rootGroup);
        }
      });
    }
  }

  onFrameAttached(dataFrame: DG.DataFrame): void {
    if (this.molColumnIdx >= 0)
      this.MoleculeColumn = this.molColumns[this.molColumnIdx].name;

    const thisViewer = this;
    this.tree.root.onscroll = async (e) => {
      await updateVisibleNodesHits(thisViewer);
    }

    this.tree.onChildNodeExpandedChanged.subscribe((group: TreeViewGroup) => updateVisibleNodesHits(thisViewer));

    this.tree.onNodeContextMenu.subscribe((args: any) => {
      const menu: DG.Menu = args.args.menu;
      const node: TreeViewGroup = args.args.item;
      menu
        .item("Add New...", () => thisViewer.openAddSketcher(node))
        .item("Edit...", () => this.openEditSketcher(node))
        .item("Remove", () => {
          node.remove();
          thisViewer.updateUI();
          thisViewer.clearFilters();
          thisViewer.TreeEncodeUpdateInProgress = true;
          thisViewer.TreeEncode = JSON.stringify(ScaffoldTreeViewer.serializeTrees(thisViewer.tree));
          thisViewer.TreeEncodeUpdateInProgress = false;
        })
    });

    this.tree.onSelectedNodeChanged.subscribe(async(node: DG.TreeViewNode) => {
       if (node.value !== null) {
         if((node.value as any).bitset === undefined)
           return;

        thisViewer.checkBoxesUpdateInProgress = true;
        this.selectGroup(node);
        thisViewer.checkBoxesUpdateInProgress = false;
        thisViewer.updateFilters();
      }
      //update the sketcher if open
      if (thisViewer.wrapper === null)
        return;

      thisViewer.wrapper.node = (node as DG.TreeViewGroup);
    });

    this.tree.onChildNodeExpandedChanged.subscribe((group: DG.TreeViewGroup) => {
      const isFolder = (group.value as any).orphans;
      if(isFolder) {
        const c = group.root.getElementsByClassName(group.expanded ? 'grok-icon fa-folder fas icon-fill' : 'grok-icon fa-folder-open fas icon-fill');
        if(c.length > 0) {
          let cl = null;
          for(let n=0; n<c.length && c[n].classList !== undefined; ++n) {
            cl = c[n].classList;
            cl.remove(group.expanded ? 'fa-folder' : 'fa-folder-open');
            cl.add(group.expanded ? 'fa-folder-open' : 'fa-folder');
          }
        }
      }
    });
    this.subs.push(dataFrame.onRowsFiltering.subscribe(() => {
      if (thisViewer.bitset != null)
        dataFrame.filter.and(thisViewer.bitset);
    }));

    this.render();
  }

  detach(): void {
    if (this.wrapper !== null) {
      this.wrapper.close()
      this.wrapper = null;
    }
    this.clearFilters();
    super.detach();
  }

  selectGroup(group: TreeViewNode) : void {
    const items = this.tree.items;
    for (let n = 0; n < items.length; ++n) {
      items[n].checked = false;
    }
    group.checked = true;
  }

  updateSizes() {
    const nodes = this.root.getElementsByClassName('d4-tree-view-node');
    for (let n = 0; n < nodes.length; ++n) {
      (nodes[n] as HTMLElement).style.height = CELL_HEIGHT + 'px';
      (nodes[n] as HTMLElement).style.width = '295px';
    }
  }

  updateUI() {
    if (this.molColumn === null) {
      this._generateLink!.style.pointerEvents = 'none';
      this._generateLink!.style.color = 'lightgrey';
      this.message = NO_MOL_COL_ERROR_MSG;
      (this._iconAdd! as any).inert = true;
      this._iconAdd!.style.color = 'grey';
      return;
    }

    const itemCount = this.tree.items.length;
    this._iconDelete!.style.visibility = itemCount > 0 ? 'visible' : 'hidden';
    this._generateLink!.style.visibility = itemCount > 0 ? 'hidden' : 'visible';
    this._message!.style.visibility = itemCount > 0 ? 'hidden' : 'visible';

    const c = this.root.getElementsByClassName('grok-icon fal fa-filter grok-icon-filter');
    if (c.length > 0) {
      (c[0] as HTMLElement).style.visibility = this.bitset === null ? 'hidden' : 'visible';
    }
  }

  render() {
    const thisViewer = this;
    this._iconAdd = ui.iconFA('plus', () => {
      thisViewer.openAddSketcher(thisViewer.tree);
    }, 'Add New Root Structure');
    this._iconAdd.style.color = "#2083d5";
    this._iconAdd.style.textAlign = "left";
    this._iconAdd.style.position = "absolute";
    this._iconAdd.style.marginLeft = "15px";
    this._iconAdd.style.fontSize = "20px";
    this.root.appendChild(this._iconAdd);

    this._iconDelete = ui.iconFA('trash-alt', () => {
      thisViewer.clear();
    }, 'Drop All Trees');
    this._iconDelete.style.color = "#2083d5";
    this._iconDelete.style.textAlign = "left";
    this._iconDelete.style.position = "absolute";
    this._iconDelete.style.marginLeft = "45px";
    this._iconDelete.style.fontSize = "19px";
    this._iconDelete.style.visibility = "hidden";
    this.root.appendChild(this._iconDelete);

    const iconClearFilter = ui.iconFA('filter', () => {
      thisViewer.clearFilters();
    }, "Clear Filter");
    iconClearFilter.classList.add('grok-icon-filter');
    iconClearFilter.style.color = "red";
    iconClearFilter.style.textAlign = "left";
    iconClearFilter.style.position = "absolute";
    iconClearFilter.style.marginLeft = "80px";
    iconClearFilter.style.fontSize = "17px";
    iconClearFilter.style.visibility = "hidden";
    this.root.appendChild(iconClearFilter);

    this.tree.root.style.position = "absolute";
    this.tree.root.style.top = "25px";
    this.tree.root.style.width = '100%';
    this.tree.root.style.height = '100%';
    this.root.appendChild(this.tree.root);

    this._message = ui.divText('', 'chem-scaffold-tree-generate-message-hint');
    this.root.appendChild(this._message);

    this._generateLink = ui.link('Generate',
      async() => await thisViewer.generateTree(),
      'Generates Scaffold Tree',
      'chem-scaffold-tree-generate-hint');
    this.root.appendChild(this._generateLink);
    this.updateSizes();
    this.updateUI();
  }

  /*
    static findOkButton(dialogRoot: HTMLElement, name: string): HTMLElement | null {
      const buttons = dialogRoot.getElementsByClassName('ui-btn ui-btn-ok');
      for (let n = 0; n < buttons.length; ++n) {
        let spans = buttons[n].getElementsByTagName("span");
        if (spans.length > 0 && spans[0].innerHTML == name)
          return buttons[n] as HTMLElement;
      }
      return null;
    }*/

  static validateNodes(childSmiles: string, parentSmiles: string): boolean {
    const parentMol = _rdKitModule.get_mol(parentSmiles);
    const parentCld = _rdKitModule.get_mol(childSmiles);

    const match: string = parentCld.get_substruct_match(parentMol);
    parentMol.delete();
    parentCld.delete();
    return match.length > 2;
  }


  static serializeTrees(treeRoot: TreeViewGroup): any {
    const json: Array<any> = [];
    for (let n = 0; n < treeRoot.children.length; ++n) {
      json[n] = ScaffoldTreeViewer.serializeTree(treeRoot.children[n] as TreeViewGroup);
    }
    console.log(json);
    return json;
  }

  static serializeTree(rootGroup: TreeViewGroup): any {
    const jsonNode: any = {};
    jsonNode.scaffold = (rootGroup.value as any).smiles;
    jsonNode.child_nodes = new Array(rootGroup.children.length);

    for (let n = 0; n < rootGroup.children.length; ++n)
      jsonNode.child_nodes[n] = ScaffoldTreeViewer.serializeTree(rootGroup.children[n] as TreeViewGroup);

    return jsonNode;
  }


  static deserializeTrees(json: any, treeRoot: TreeViewGroup, createGroup: Function) : number {
    let countNodes = 0;
    for (let n = 0; n < json.length; ++n) {
      countNodes += ScaffoldTreeViewer.deserializeTree(json[n], treeRoot, (molStr: string, rootGroup: TreeViewGroup,  countNodes: number) => {
        return createGroup(molStr, rootGroup, countNodes);
      }, 0);
    }
    return countNodes;
  }

  static deserializeTree(json: any, rootGroup: TreeViewGroup, createGroup: Function, countNodes: number) : number {
    const molStr = json.scaffold;
    if (molStr === null || molStr === undefined) {
      console.error('Scaffold is null or undefined.')
      return countNodes;
    }

    const group: TreeViewGroup = createGroup(molStr, rootGroup, countNodes);
    if (group === null)
      return countNodes;

    if (json.child_nodes === undefined)
      json.child_nodes = [];

    //console.log('child_nodes: ' + json.child_nodes.length);
    for (let n = 0; n < json.child_nodes.length; ++n) {
      //console.log('processing child: ' + n);
      countNodes += ScaffoldTreeViewer.deserializeTree(json.child_nodes[n], group, createGroup, countNodes);
    }

    ++countNodes;

    if (group.children.length === 0)
      enableNodeExtendArrow(group, false);

    return countNodes;
  }
}

class SketcherDialogWrapper {
  private readonly dialog: DG.Dialog;
  private readonly sketcher: Sketcher;
  private success: boolean;
  private group: DG.TreeViewGroup;
  private isMolBlock: boolean;
  private activeElement : HTMLElement | null = null;

  constructor(title: string, actionName: string, group: DG.TreeViewGroup, action: Function, validate: Function, onCancel: Function | null = null, onClose: Function | null = null, onStrucChanged: Function | null = null) {
    this.success = true;
    this.dialog = ui.dialog({title: title});
    this.group = group;
    const v = this.group.value as any;
    const molStr = v === null ? '' : v.smiles;
    this.isMolBlock = isMolBlock(molStr);

    const thisWrapper = this;

    const validLabel = ui.label(SketcherDialogWrapper.validationMessage(true, true));
    validLabel.style.height = "30px";
    validLabel.style.color = SketcherDialogWrapper.validationColor(true);

    this.sketcher = new Sketcher();
    this.isMolBlock ? this.sketcher.setMolFile(molStr) : this.sketcher.setSmiles(molStr);

    const molStrTmp = thisWrapper.sketcher.getMolFile();
    let valid = validate(molStrTmp, thisWrapper.node);
    validLabel.style.color = SketcherDialogWrapper.validationColor(valid);
    validLabel.innerText = SketcherDialogWrapper.validationMessage(valid, true);

    this.sketcher.onChanged.subscribe(() => {
      const molStr = thisWrapper.sketcher.getMolFile();//thisWrapper.isMolBlock ? thisWrapper.sketcher.getMolFile() : thisWrapper.sketcher.getSmiles();
      valid = validate(molStr, thisWrapper.node);
      validLabel.style.color = SketcherDialogWrapper.validationColor(valid);
      validLabel.innerText = SketcherDialogWrapper.validationMessage(valid, true);

      if (onStrucChanged !== null)
        onStrucChanged(molStr, thisWrapper.sketcher.getMolFile());
      //const btnOk : HTMLElement | null = thisWrapper.dialog === null ? null : ScaffoldTreeViewer.findOkButton(thisWrapper.dialog.root, actionName);
      //(btnOk as any).disabled = !success;
    });

    this.dialog.add(this.sketcher);
    this.dialog.add(validLabel);
    this.dialog.addButton("Reset", () => {
      thisWrapper.isMolBlock ? thisWrapper.sketcher.setMolFile(molStr) : thisWrapper.sketcher.setSmiles(molStr);
    });
    this.dialog.addButton(actionName, () => {
      const molStr = thisWrapper.sketcher.getMolFile();//thisWrapper.isMolBlock ? thisWrapper.sketcher.getMolFile() : thisWrapper.sketcher.getSmiles();
      action(molStr, thisWrapper.node, valid);
    });

    if (onCancel != null)
      this.dialog.onCancel(onCancel);

    if (onClose != null)
      this.dialog.onCancel(onClose);
  }

  set node(node: DG.TreeViewGroup) {
    this.group = node;
    const v = node.value as any;
    const molStr = v === null ? '' : v.smiles;
    this.isMolBlock = isMolBlock(molStr);
    this.isMolBlock ? this.sketcher.setMolFile(molStr) : this.sketcher.setSmiles(molStr);
  }

  get node() {
    return this.group;
  }

  show(): void {
     this.activeElement = document.activeElement instanceof HTMLElement ? document.activeElement as HTMLElement : null;
     this.dialog?.show();
  }

  close(): void {
    this.activeElement!.style!.display = 'none';
    this.dialog?.close();
    const thisWrapper = this;
    setTimeout(() => thisWrapper.activeElement!.style!.display = '',20);
  }

  static create(title: string, actionName: string, group: DG.TreeViewGroup, action: Function, validate: Function, onCancel: Function | null = null, onClose: Function | null = null, onStrucChanged: Function | null = null): SketcherDialogWrapper {
    return new SketcherDialogWrapper(title, actionName, group, action, validate, onCancel, onClose, onStrucChanged);
  }

  static validationMessage(success: boolean, parentCheck: boolean): string {
    if (success)
      return '';

    return parentCheck ? "The edited molecule is not a superstructure of its parent" :
      "The edited molecule is not a substructure of its children";
  }

  static validationColor(success: boolean): string {
    return success ? SketcherDialogWrapper.SUCCESS_MSG_COLOR : SketcherDialogWrapper.FAILURE_MSG_COLOR;
  }

  static SUCCESS_MSG_COLOR = 'green';
  static FAILURE_MSG_COLOR = 'red';
}