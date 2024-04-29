import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Template, TEMPLATES_FOLDER } from './constants';
import { createConditionalInput, createInputForProperty, createLinearInput } from './admetox-utils';

export class AdmeticaBaseEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
	templatesInput: DG.InputBase; 
  addPiechartInput = ui.boolInput('Add piechart', true);
  modelsSettingsDiv = ui.inputs([]);
  modelsSettingsIcon: HTMLElement;

  constructor() {
    this.tableInput = ui.tableInput('Table', grok.shell.tv.dataFrame, grok.shell.tables, () => {
      // TODO: add onTableAttached
    });
		this.colInput = ui.columnInput('Column', grok.shell.tv.dataFrame, grok.shell.tv.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE));
		this.templatesInput = ui.choiceInput('Templates', 'template', ['template'], async () =>  {
			if(settingsOpened)
				await this.createModelsSettingsDiv(this.modelsSettingsDiv);
		});
    this.modelsSettingsIcon = ui.icons.settings(async () => {
			settingsOpened = !settingsOpened;
      if (settingsOpened)
				await this.createModelsSettingsDiv(this.modelsSettingsDiv);
      else
			  ui.empty(this.modelsSettingsDiv);  
    }, 'Modify template parameters');
		this.templatesInput.root.prepend(this.modelsSettingsIcon);
		let settingsOpened = false;
  }

	private async getTemplates(): Promise<string[]> {
    const files = await grok.dapi.files.list(TEMPLATES_FOLDER);
    return files.map((file) => file.fileName.split('.')[0]);
	}

  private async createModelsSettingsDiv(paramsForm: HTMLElement): Promise<HTMLElement> {
		const items = await grok.dapi.files.list(TEMPLATES_FOLDER);
		const fileName = items[0].fileName;
		const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
		const properties = JSON.parse(propertiesJson);
		ui.empty(paramsForm);
		const treeControl = this.createTreeControl(properties);
		paramsForm.appendChild(treeControl);
		return paramsForm;
	}
	

	private createTreeControl(template: Template) {
		const treeView = ui.tree();
		treeView.root.style.overflow = 'hidden';
		const inputs = ui.divV([]);

		const createInputs = (model: any) => {
			if (!model) {
				ui.empty(inputs);
				return;
			}

			ui.empty(inputs);
			inputs.classList.add('admetox-input-form');
			const properties = model.properties;
			const coloring = model.coloring;
	
			properties.forEach((p: any) => {
				const input = createInputForProperty(p);
				if (input)
					inputs.appendChild(input);
			});
	
			if (coloring.type === 'Conditional') {
				const conditionalInput = createConditionalInput(coloring, model);
				inputs.appendChild(conditionalInput);
			} else if (coloring.type === 'Linear') {
				const linearInput = createLinearInput(coloring);
				inputs.appendChild(linearInput);
			}
		}
		template.subgroup.forEach((subgroup: any) => {
			const groupNode = treeView.group(subgroup.name);
			groupNode.expanded = false;
			groupNode.enableCheckBox(false);
	
			subgroup.models.forEach((model: any) => {
				const modelNode = groupNode.item(model.name);
				modelNode.checked = false;
				modelNode.enableCheckBox(false);
			});
		});

		treeView.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
			const nodeName = node.text;
			const model = template.subgroup
			  .flatMap(subgroup => subgroup.models)
				.find(model => nodeName.includes(model.name));
			createInputs(model);
		})

		treeView.onNodeMouseEnter.subscribe((node: DG.TreeViewNode) => {
			const nodeName = node.text;
			const model = template.subgroup
			  .flatMap(subgroup => subgroup.models)
				.find(model => nodeName.includes(model.name));
			if (!model) return;
			node.root.onmouseenter = (e) => ui.tooltip.show(model.units, e.x, e.y);
			node.root.onmouseleave = (e) => ui.tooltip.hide();
		});

		const tabContent = ui.divH([treeView.root, inputs]);
		tabContent.style.gap = '10px';
		return tabContent;
	}
	

  public getEditor(): HTMLElement {
    return ui.div([
      this.tableInput.root,
      this.colInput ? this.colInput.root : null,
			this.templatesInput,
			this.modelsSettingsDiv,
      this.addPiechartInput.root,
    ], { style: { minWidth: '450px' }, classes: 'ui-form' });
  }

  public getParams() {
    return {
      table: this.tableInput.value!,
      col: this.colInput ? this.colInput.value! : null,
			templatesName: this.templatesInput.value,
    };
  }
}
