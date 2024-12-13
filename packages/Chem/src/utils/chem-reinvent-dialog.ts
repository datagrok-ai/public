import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getConfigFiles } from '../package';

export class ReinventBaseEditor {
  molInput: DG.InputBase<string | null>;
  targetInput: DG.ChoiceInput<string>;
  tabControl: DG.TabControl;
  tree: DG.TreeViewGroup = ui.tree();

  constructor() {
    this.molInput = ui.input.molecule('Ligand', { value: "OC(CN1CCCC1)NC(CCC1)CC1Cl" });
    this.targetInput = ui.input.choice('Target', { nullable: false }) as DG.ChoiceInput<string>;
    this.tabControl = ui.tabControl(null, true);
    this.tabControl.addPane('Binding affinity', () => ui.panel([this.targetInput.root]));
    this.tabControl.addPane('ADME', () => {
        const container = ui.div();
        this.admeticaTree().then((admeTree) => {
          container.appendChild(ui.panel([admeTree.root]));
        });
        return container;
    });
    this.initTargets();
  }

  private async admeticaTree() {
    this.tree.root.style.overflow = 'hidden';
    const properties = ['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity'];
    for (const prop of properties) {
        const group = this.tree.group(prop);
        group.enableCheckBox(false);
        group.expanded = true;

        const models = await grok.functions.call('Admetica:getModels', {property: prop});
        for (const model of models) {
            const item = group.item(model);
            item.enableCheckBox(false);
        }
    }
    return this.tree;
  }

  private async initTargets() {
    const targets = await getConfigFiles();
    this.targetInput.items = targets;
  }

  public getEditor(): HTMLElement {
    const moleculeBlock = ui.panel([
        ui.h1('Molecule input'),
        this.molInput.root
    ]);

    const parametersBlock = ui.panel([
        ui.h1('Parameter settings'),
        this.tabControl.root
    ]);

    return ui.divV([moleculeBlock, parametersBlock], {
      style: { minWidth: '450px', display: 'flex', flexDirection: 'column', gap: '15px' },
      classes: 'ui-form',
    });
  }

  public getParams() {
    return {
      ligand: this.molInput.value!,
      target: this.targetInput.value!,
      models: this.tree?.items.filter((item) => !(item instanceof DG.TreeViewGroup) && item.checked).map((item) => item.text),
    };
  }
}