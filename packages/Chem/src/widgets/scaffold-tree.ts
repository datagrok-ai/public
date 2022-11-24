import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from "cash-dom";
import {_rdKitModule } from '../utils/chem-common-rdkit';
import {chem} from "datagrok-api/grok";
import {TreeViewGroup} from "datagrok-api/dg";
import Sketcher = chem.Sketcher;
import {RDKitCellRenderer} from "../rendering/rdkit-cell-renderer";

function renderMolecule(molStr: string, width: number, height: number): HTMLElement {
    const moleculeHost = ui.canvas(width, height);
    $(moleculeHost).addClass('chem-canvas');
    const r = window.devicePixelRatio;
    moleculeHost.width = width * r;
    moleculeHost.height = height * r;
    moleculeHost.style.width = width.toString() + 'px';
    moleculeHost.style.height = height.toString() + 'px';
    const renderer = new RDKitCellRenderer(_rdKitModule);

    // @ts-ignore
    renderer.render(moleculeHost.getContext('2d')!, 0, 0, width, height, DG.GridCell.fromValue(molStr));
    return ui.divV([moleculeHost], 'chem-mol-box');
}

function createMolDiv(smiles : string) : HTMLElement{
    return renderMolecule(smiles, 75,75);
}

export class ScaffoldTreeViewer extends DG.JsViewer {
    private tree : DG.TreeViewGroup;

    constructor() {
      super();
      this.tree = ui.tree();
    }

    onFrameAttached(dataFrame: DG.DataFrame): void {
      const thisViewer = this;
      let wrapper : SketcherDialogWrapper | null = null;
      let h = this.tree.onNodeContextMenu.subscribe((args : any) => {
        const menu : DG.Menu = args.args.menu;
        const node : TreeViewGroup = args.args.item;
        const smiles = (node.value as any).smiles;
        menu.item("Add New...", () => {
        wrapper = SketcherDialogWrapper.create("Add New Scaffold...", "Add", node, (smilesSketcher: string, node : TreeViewGroup) => {
           node.group(createMolDiv(smilesSketcher), {smiles : smilesSketcher});
           thisViewer.updateSizes();
           wrapper?.close();
           wrapper = null;
          }, (smilesSketcher: string) => {
           return ScaffoldTreeViewer.validateNodes(smilesSketcher, smiles);
         });
          wrapper.show();
        });

        menu.item("Edit...", () => {
          wrapper = SketcherDialogWrapper.create("Edit Scaffold...", "Save", node, (smilesSketcher: string, node : TreeViewGroup) => {
           while (node.captionLabel.firstChild) {
            node.captionLabel.removeChild(node.captionLabel.firstChild);
           }
           node.captionLabel.appendChild(createMolDiv(smilesSketcher));
           node.value = {smiles : smilesSketcher};
           wrapper?.close();
           wrapper = null;
          }, (smilesSketcher: string, node : TreeViewGroup) => {

            if(node.parent === null)
              return true;

            let success = true;
            if(node.parent.value !== null) {
              const smilesParent = (node.parent.value as any).smiles;
              success = ScaffoldTreeViewer.validateNodes(smilesSketcher, smilesParent);
              if (!success)
               return false;
            }
            const children = node.items;
            for(let n=0; n<children.length; ++n) {
              success = ScaffoldTreeViewer.validateNodes((children[n].value as any).smiles, smilesSketcher);
              if(!success)
               return false;
            }

            return true;
          });

          wrapper.show();
        });

        menu.item("Remove", () => {
          node.remove();
        });
      });
      this.subs.push(h);

      h = this.tree.onSelectedNodeChanged.subscribe((node : DG.TreeViewNode) => {
       if(wrapper === null)
        return;

       wrapper.node = (node as DG.TreeViewGroup);
      });
      this.subs.push(h);

      const scaffold = "O=C1CN=C(c2ccccc2N1)C3CCCCC3";
      this.tree.group(createMolDiv(scaffold), {smiles : scaffold});
      this.render();
    }

    detach(): void {
      super.detach();
    }

    updateSizes() {
      let nodes = this.root.getElementsByClassName('d4-tree-view-node');
      for (let n=0; n<nodes.length; ++n) {
        (nodes[n] as HTMLElement).style.height = "75px";
      }
    }

    render() {
      this.root.appendChild(this.tree.root);
      this.updateSizes();
    }

    static findOkButton(dialogRoot : HTMLElement, name : string) : HTMLElement | null {
      const buttons = dialogRoot.getElementsByClassName('ui-btn ui-btn-ok');
      for (let n=0; n<buttons.length; ++n) {
        let spans = buttons[n].getElementsByTagName("span");
        if(spans.length > 0 && spans[0].innerHTML == name)
         return buttons[n] as HTMLElement;
     }
     return null;
    }

    static validateNodes(childSmiles : string, parentSmiles: string) : boolean {
      const parentMol = _rdKitModule.get_mol(parentSmiles);
      const parentCld = _rdKitModule.get_mol(childSmiles);
      const match : string = parentCld.get_substruct_match(parentMol);
      parentMol.delete();
      parentCld.delete();
      return match.length > 2;
    }
}

class SketcherDialogWrapper {
  private readonly dialog : DG.Dialog;
  private readonly sketcher: Sketcher;
  private success: boolean;
  private group : DG.TreeViewGroup;

  constructor(title : string, actionName: string, group : DG.TreeViewGroup, action : Function, validate : Function) {
    this.success = true;
    this.dialog = ui.dialog({title : title});//  dialog = ui.dialog({title: "Add New Scaffold..."});
    this.group = group;
    const smiles = (this.group.value as any).smiles;

    const validLabel = ui.label("The Node Structure is Valid.");
    validLabel.style.height = "30px";
    validLabel.style.color = "green";
    this.dialog.add(validLabel);

    this.sketcher = new Sketcher();
    this.dialog.add(ui.button("Reset Structure", () => {
     this.sketcher?.setSmiles(smiles);
    }));

    this.sketcher.setSmiles(smiles);
    this.sketcher.onChanged.subscribe(() => {
      const success = validate(this.sketcher.getSmiles(), this.node);
      validLabel.style.color = success ? "green" : "red";
      validLabel.innerText =  success ? "The Node Structure is Valid." : "The Node Structure is Invalid.";
      const btnOk : HTMLElement | null = this.dialog === null ? null : ScaffoldTreeViewer.findOkButton(this.dialog.root, actionName);
      (btnOk as any).disabled = !success;
    });

    this.dialog.add(this.sketcher);
    this.dialog.addButton(actionName, () => {
     action(this.sketcher.getSmiles(), this.node);
    });
  }

  set node(node : DG.TreeViewGroup) {
    this.group = node;
    this.sketcher.setSmiles((node.value as any).smiles);
  }

  get node() {return this.group;}
  show() : void {this.dialog?.show();}
  close() : void {this.dialog?.close();}

  static create(title : string, actionName: string, group : DG.TreeViewGroup, action : Function, validate : Function) : SketcherDialogWrapper {
    return new SketcherDialogWrapper(title, actionName, group, action, validate);
  }
}