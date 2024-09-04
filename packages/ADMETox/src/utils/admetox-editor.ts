import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { Model, ModelColoring, Subgroup, Template, TEMPLATES_FOLDER } from './constants';

import '../css/admetox.css';

export class AdmeticaBaseEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot: HTMLElement;
  templatesInput: DG.ChoiceInput<string>;
  saveButton: HTMLElement;
  addPiechartInput = ui.input.bool('Add piechart', {value: true});
  modelsSettingsDiv = ui.inputs([]);
  modelsSettingsIcon: HTMLElement;
  expanded: boolean = false;
  tree: DG.TreeViewGroup = ui.tree();
  properties: Template | null = null;
  updatedProperties: Template | null = null;

  constructor() {
    this.tableInput = ui.input.table('Table', {value: grok.shell.tv.dataFrame, onValueChanged: () => this.onTableInputChanged()});
    this.colInput = ui.input.column('Molecules', {
      table: grok.shell.tv.dataFrame,
      filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE,
      value: grok.shell.tv.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!
    });
    this.colInputRoot = this.colInput.root;
    this.saveButton = ui.button('Save...', () => {
      this.dialogForTemplateSave();
    });
    this.saveButton.classList.add('admetox-save-button');
    
    this.templatesInput = ui.input.choice('Template', {onValueChanged: async () =>  {
      if (settingsOpened)
        await this.createModelsSettingsDiv(this.modelsSettingsDiv);
    }}) as DG.ChoiceInput<string>;

    this.modelsSettingsIcon = ui.icons.settings(async () => {
      settingsOpened = !settingsOpened;
      if (settingsOpened)
        await this.createModelsSettingsDiv(this.modelsSettingsDiv);
      else
        ui.empty(this.modelsSettingsDiv);
    }, 'Modify template parameters');
    this.templatesInput.root.prepend(this.modelsSettingsIcon);
    let settingsOpened = false;
    this.initTemplates();
  }
  
  private async initTemplates() {
    const templates = await this.getTemplates();
    this.templatesInput.items = templates;
    this.templatesInput.value = templates[0];
  }
  
  private async getTemplates(): Promise<string[]> {
    const files = await grok.dapi.files.list(TEMPLATES_FOLDER);
    return files.map((file) => file.fileName.split('.')[0]);
  }
  
  private async getPropertiesFromTemplate() {
    const fileName = `${this.templatesInput.value}.json`;
    const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
    return JSON.parse(propertiesJson);
  }

  private async dialogForTemplateSave(): Promise<void> {
    const templateNameInput = ui.input.string('Name');
    ui.dialog({ title: 'Save template' })
      .add(templateNameInput)
      .onOK(async () => {
        const templateName = templateNameInput.value;
        await grok.dapi.files.writeAsText(`${TEMPLATES_FOLDER}/${templateName}.json`, JSON.stringify(this.updatedProperties));
        grok.shell.info(`Template "${templateName}" has been successfully saved!`);
        //await this.initTemplates();
      })
      .show();
  }
  
  private async createModelsSettingsDiv(paramsForm: HTMLElement): Promise<HTMLElement> {
    this.properties = await this.getPropertiesFromTemplate();
    this.updatedProperties = JSON.parse(JSON.stringify(this.properties));
    ui.empty(paramsForm);
    this.tree = ui.tree();
    this.expanded = false;
    const treeControl = this.createTreeControl(this.updatedProperties!);
    paramsForm.appendChild(treeControl);
    this.tree.currentItem = this.tree.items.find(item => item instanceof DG.TreeViewNode)!;
    return paramsForm;
  }

  private createInputForProperty(property: any) {
    if (property.property.skip) {
     return; 
    }
    const object = property.property.inputType === DG.InputType.Map ? {} : property.object;
    const prop = DG.Property.fromOptions(property.property);
    const input = DG.InputBase.forProperty(prop, object);
    if (!property.property.enable) {
      input.root.style.pointerEvents = 'none';
    }
    const key = property.property.name as keyof typeof property.object;
    input.value = property.object[key];
    input.addCaption('');
    input.onChanged.subscribe(() => {
      property.object[key] = input.value;
      const areEqual = JSON.stringify(this.properties) === JSON.stringify(this.updatedProperties);
      this.saveButton.style.visibility = areEqual ? 'hidden' : 'visible';
    });
    return input.root;
  }
  
  private createConditionalInput(coloring: ModelColoring, model: Model) {
    const conditionalColors = Object.entries(coloring).slice(1);
    const conditionalColoring = `{${conditionalColors.map(([range, color]) => `"${range}":"${color}"`).join(",")}}`;
    const patternsInp = ui.patternsInput(JSON.parse(conditionalColoring));
    const inputs = patternsInp.querySelectorAll('.ui-input-editor');
    const rangesProperty = model.properties.find((prop: any) => prop.property.name === 'ranges');
    inputs.forEach((input, index) => input.addEventListener('mousemove', function(event) {
      const mouseEvent = event as MouseEvent;
      const key = conditionalColors[index][0];
      let value = '';
      if (rangesProperty)
        value = rangesProperty.object.ranges[key];
      ui.tooltip.show(value, mouseEvent.x, mouseEvent.y);
    }));
    return patternsInp;
  }
  
  private createLinearInput(coloring: ModelColoring) {
    const linearInput = ui.schemeInput(JSON.parse(coloring.colors!) as number[]);
    linearInput.removeChild(linearInput.firstChild!);
    const div = ui.divV([linearInput]);
    linearInput.style.pointerEvents = 'none';
    return div;
  }
  
  private createInputsForCategories(group: Subgroup, inputs: HTMLElement): void{
    ui.empty(inputs);
    inputs.classList.add('admetox-input-form');
    inputs.appendChild(ui.divText(group.name, 'admetox-descriptor-name'));
    inputs.appendChild(ui.input.textArea('', {size: {width: 100, height: 150}, value: group.description}).root);
  }
  
  private createInputsForModels(model: Model, inputs: HTMLElement): void {
    ui.empty(inputs);
    inputs.classList.add('admetox-input-form');
    inputs.appendChild(ui.divText(model.name, 'admetox-descriptor-name'));
    model.properties.forEach(p => {
      const input = this.createInputForProperty(p);
      if (input) inputs.appendChild(input);
    });
    
    const coloring = model.coloring;
    if (coloring.type === 'Conditional')
      inputs.appendChild(this.createConditionalInput(coloring, model));
    else if (coloring.type === DG.COLOR_CODING_TYPE.LINEAR)
      inputs.appendChild(this.createLinearInput(coloring));
  }
  
  private createTreeGroup(template: Template): void {
    template.subgroup.forEach((subgroup: Subgroup) => {
      const groupNode = this.tree.group(subgroup.name);
      groupNode.expanded = !this.expanded;
      this.expanded = true;
      groupNode.enableCheckBox(false);
      
      subgroup.models.forEach((model: Model) => {
        const modelNode = groupNode.item(model.name);
        modelNode.checked = false;
        modelNode.enableCheckBox(false);
      });
    });
  }
  
  private getModel(template: Template, nodeName: string): Model | null {
    const model = template.subgroup
      .flatMap(subgroup => subgroup.models)
      .find(model => nodeName === model.name);
    if (model) return model;
    return null;
  }
  
  private createTreeControl(template: Template): HTMLDivElement {
    this.createTreeGroup(template);
    this.tree.root.classList.add('admetox-tree');
    const inputs = ui.divV([]);
    
    this.tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
      const subgroup = template.subgroup.find(subgroup => subgroup.name === node.text);
      const model = this.getModel(template, node.text)
      if (subgroup)
        this.createInputsForCategories(subgroup, inputs);
      else if (model)
        this.createInputsForModels(model, inputs);
    });
    
    this.tree.onNodeMouseEnter.subscribe((node: DG.TreeViewNode) => {
      const model = this.getModel(template, node.text);
      if (!model) return;
      node.root.onmouseenter = (e) => ui.tooltip.show(model.units, e.x, e.y);
      node.root.onmouseleave = (e) => ui.tooltip.hide();
    });
    
    const tabContent = ui.divH([this.tree.root, inputs]);
    tabContent.classList.add('admetox-tab-content');
    return tabContent;
  }
  
  public getEditor(): HTMLElement {
    return ui.div([
      this.tableInput,
      this.colInputRoot,
      ui.divH([this.templatesInput.root, this.saveButton]),
      this.modelsSettingsDiv,
      this.addPiechartInput.root,
    ], { style: { minWidth: '450px' }, classes: 'ui-form' });
  }
  
  public getParams() {
    return {
      table: this.tableInput.value!,
      col: this.colInput ? this.colInput.value! : null,
      templatesName: this.templatesInput.value,
      models: this.tree?.items.filter((item) => !(item instanceof DG.TreeViewGroup) && item.checked).map((item) => item.text),
      addPiechart: this.addPiechartInput.value
    };
  }
  
  onTableInputChanged(): void {
    this.colInput = ui.input.column('Molecule', {
      table: this.tableInput.value!,
      value: this.tableInput.value!.columns.bySemType(DG.SEMTYPE.MOLECULE)!,
      filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE
    });
    ui.empty(this.colInputRoot);
    Array.from(this.colInput.root.children).forEach((it) => this.colInputRoot.append(it));
  }
}