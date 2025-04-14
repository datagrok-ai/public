import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { Model, ModelColoring, Subgroup, Template, TEMPLATES_FOLDER } from './constants';
import { getTemplates } from './admetica-utils';
import '../css/admetica.css';

export class AdmeticaBaseEditor {
  tableInput: DG.InputBase<DG.DataFrame | null>;
  colInput!: DG.InputBase<DG.Column | null>;
  colInputRoot: HTMLElement;
  templatesInput: DG.ChoiceInput<string>;
  folderIcon: HTMLElement;
  saveButton: HTMLElement;
  addPiechartInput = ui.input.bool('Add piechart', {value: true});
  addFormInput = ui.input.bool('Add form', {value: true});
  modelsSettingsDiv = ui.inputs([]);
  tree: DG.TreeViewGroup = ui.tree();
  properties: Template | null = null;
  updatedProperties: Template | null = null;
  currentItem: DG.TreeViewNode | null = null;

  constructor() {
    this.tableInput = ui.input.table('Table', {value: grok.shell.tv.dataFrame, onValueChanged: () => this.onTableInputChanged()});
    this.colInput = ui.input.column('Molecules', {
      table: grok.shell.tv.dataFrame,
      filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE,
      value: grok.shell.tv.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE)!,
      nullable: false
    });
    this.colInputRoot = this.colInput.root;
    this.saveButton = ui.button('Save...', () => {
      this.dialogForTemplateSave();
    });
    this.saveButton.classList.add('admetica-save-button');
    
    this.templatesInput = ui.input.choice('Template', {
      onValueChanged: async () =>  await this.createModelsSettingsDiv(this.modelsSettingsDiv),
      nullable: false
    }) as DG.ChoiceInput<string>;

    this.folderIcon = ui.iconFA('folder', () => {
      ui.dialog({title:'Manage files'})
        .add(ui.fileBrowser({path: TEMPLATES_FOLDER}).root)
        .onOK(async () => await this.initTemplates())
        .show();
    }, 'Manage templates (rename/delete)');
    this.folderIcon.style.marginLeft = '10px';
    this.templatesInput.root.append(this.folderIcon);

    this.initTemplates();
  }
  
  private async initTemplates(template?: string) {
    const templates = await getTemplates();
    this.templatesInput.items = templates;
    if (template)
      this.templatesInput.value = template;
  }
  
  private async getPropertiesFromTemplate() {
    const fileName = `${this.templatesInput.value}.json`;
    const propertiesJson = await grok.dapi.files.readAsText(`${TEMPLATES_FOLDER}/${fileName}`);
    return JSON.parse(propertiesJson);
  }

  private async dialogForTemplateSave(): Promise<void> {
    const templateNameInput = ui.input.string('Name');
    templateNameInput.addValidator(
      (template: string) => this.templatesInput.items.includes(template) ? 'Template with this name already exists' : null);
    ui.dialog({ title: 'Save template' })
      .add(templateNameInput)
      .onOK(async () => {
        const templateName = templateNameInput.value;
        await grok.dapi.files.writeAsText(`${TEMPLATES_FOLDER}/${templateName}.json`, JSON.stringify(this.updatedProperties));
        grok.shell.info(`Template "${templateName}" has been successfully saved!`);
        await this.initTemplates(templateName);
      })
      .show();
  }
  
  private async createModelsSettingsDiv(paramsForm: HTMLElement): Promise<HTMLElement> {
    this.properties = await this.getPropertiesFromTemplate();
    this.updatedProperties = JSON.parse(JSON.stringify(this.properties));
    ui.empty(paramsForm);
    this.tree = ui.tree();
    const treeControl = this.createTreeControl(this.updatedProperties!);
    paramsForm.appendChild(treeControl);
    this.tree.currentItem = this.currentItem ?? this.tree.items.find(item => item instanceof DG.TreeViewNode)!;
    return paramsForm;
  }

  private createInputForProperty(property: any, model: Model) {
    if (property.property.skip) {
     return; 
    }
    const object = property.property.inputType === DG.InputType.Map ? {} : property.object;
    const prop = DG.Property.fromOptions(property.property);
    const input = DG.InputBase.forProperty(prop, object);
    const isLinear = model.coloring.type === DG.COLOR_CODING_TYPE.LINEAR;

    if (!property.property.enable)
      input.root.style.pointerEvents = 'none';

    const key = property.property.name as keyof typeof property.object;
    input.value = property.object[key];
    
    if (input.inputType === DG.InputType.TextArea)
      input.root.style.height = `${this.getTextHeight(property.object[key])}px`;

    input.addCaption('');
    input.onChanged.subscribe(() => {
      property.object[key] = input.value;
      if (key === 'direction' && isLinear)
        [model.coloring.min, model.coloring.max] = [model.coloring.max, model.coloring.min];  
      const areEqual = JSON.stringify(this.properties) === JSON.stringify(this.updatedProperties);
      this.saveButton.style.visibility = areEqual ? 'hidden' : 'visible';
    });

    if (input.inputType === DG.InputType.Map) {
      const table = input.root.querySelector('.d4-item-table') as HTMLTableElement;
      Array.from(table!.rows).forEach(row => row.lastElementChild?.remove());
      table.style.pointerEvents = 'none';
      return table;
    }

    return input.root;
  }

  private createConditionalInput(coloring: ModelColoring, model: Model) {
    let conditionalColors = Object.entries(coloring).slice(1);
    const conditionalColoring = `{${conditionalColors.map(([range, color]) => `"${range}":"${color}"`).join(",")}}`;
    const patternsInp = ui.patternsInput(JSON.parse(conditionalColoring));
    const inputs = patternsInp.querySelectorAll('.ui-input-editor');
    const colorBars = patternsInp.querySelectorAll('.d4-color-bar');
    
    let rangesProperty = model.properties.find((prop: any) => prop.property.name === 'ranges');
    let coloringProperty = model.coloring as any;
    
    inputs.forEach((input, index) => {
      let key = conditionalColors[index][0];
      const inputElement = input as HTMLInputElement;
      
      inputElement.addEventListener('mousemove', (event) => {
        const mouseEvent = event as MouseEvent;
        const value = rangesProperty?.object.ranges[key] || '';
        ui.tooltip.show(value, mouseEvent.x, mouseEvent.y);
      });
  
      inputElement.addEventListener('input', () => {
        const newValue = inputElement.value;
        
        if (rangesProperty) {
          const ranges = rangesProperty.object.ranges;
          const oldValue = ranges[key];
          delete ranges[key];
          ranges[newValue] = oldValue;
        }
  
        if (coloringProperty) {
          const oldColorValue = coloringProperty[key];
          delete coloringProperty[key];
          coloringProperty[newValue] = oldColorValue;
        }
  
        conditionalColors[index][0] = newValue;
        key = newValue;

        const areEqual = JSON.stringify(this.properties) === JSON.stringify(this.updatedProperties);
        this.saveButton.style.visibility = areEqual ? 'hidden' : 'visible';
      });
    });

    colorBars.forEach((colorBar, index) => {
  
      const handleColorChange = () => {
        let key = conditionalColors[index][0];
        const newColor = window.getComputedStyle(colorBar).backgroundColor;
        const hexColor = this.rgbToHex(newColor);
  
        if (coloringProperty) {
          coloringProperty[key] = hexColor;
        }
  
        const areEqual = JSON.stringify(this.properties) === JSON.stringify(this.updatedProperties);
        this.saveButton.style.visibility = areEqual ? 'hidden' : 'visible';
      };
  
      const observer = new MutationObserver(() => {
        handleColorChange();
      });
  
      observer.observe(colorBar, { attributes: true, attributeFilter: ['style'] });
    });
  
    return patternsInp;
  }
  
  private rgbToHex(rgb: string): string {
    const result = rgb.match(/\d+/g);
    if (!result || result.length < 3) return '#000000';
  
    const r = parseInt(result[0]).toString(16).padStart(2, '0');
    const g = parseInt(result[1]).toString(16).padStart(2, '0');
    const b = parseInt(result[2]).toString(16).padStart(2, '0');
  
    return `#${r}${g}${b}`;
  }  
  
  private createLinearInput(coloring: ModelColoring) {
    const linearInput = ui.schemeInput(JSON.parse(coloring.colors!) as number[]);
    const canvas = linearInput.querySelector('canvas');
    canvas!.style.width = '100px';
    canvas!.style.height = '20px';
    const div = ui.divV([canvas]);
    linearInput.style.pointerEvents = 'none';
    return div;
  }
  
  private createInputsForCategories(group: Subgroup, inputs: HTMLElement): void{
    ui.empty(inputs);
    inputs.classList.add('admetica-input-form');
    inputs.appendChild(ui.divText(group.name, 'admetica-descriptor-name'));
    inputs.appendChild(ui.input.textArea('', {size: {width: 100, height: 150}, value: group.description}).root);
  }
  
  private createInputsForModels(model: Model, inputs: HTMLElement): void {
    ui.empty(inputs);
    inputs.classList.add('admetica-input-form');
    inputs.appendChild(ui.divText(model.name, 'admetica-descriptor-name'));
    model.properties.forEach(p => {
      const input = this.createInputForProperty(p, model);
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
      groupNode.expanded = subgroup.expanded;
      groupNode.enableCheckBox(subgroup.checked);
      
      subgroup.models.forEach((model: Model) => {
        const modelNode = groupNode.item(model.name);
        modelNode.checked = model.checked;
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
    this.tree.root.classList.add('admetica-tree');
    const inputs = ui.divV([]);
    
    this.tree.onSelectedNodeChanged.subscribe((node: DG.TreeViewNode) => {
      this.currentItem = node;
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

    this.tree.onChildNodeExpandedChanged.subscribe((node: DG.TreeViewGroup) => {
      const subgroup = template.subgroup.find(subgroup => subgroup?.name === node.text);
      if (subgroup)
        subgroup.expanded = node.expanded;
    });

    this.tree.onNodeCheckBoxToggled.subscribe((node) => {
      if (node instanceof DG.TreeViewGroup) {
        const subgroup = template.subgroup.find(subgroup => subgroup?.name === node.text);
        subgroup!.checked = node.checked;
        subgroup!.models.forEach((model: Model) => {
          model.checked = node.checked;
        })
      } else {
        const model = this.getModel(template, node.text);
        model!.checked = node.checked;
      }
    });
    
    const tabContent = ui.divH([this.tree.root, inputs]);
    tabContent.classList.add('admetica-tab-content');
    return tabContent;
  }
  
  public getEditor(): HTMLElement {
    return ui.div([
      this.tableInput,
      this.colInputRoot,
      ui.divH([this.templatesInput.root, this.saveButton]),
      this.modelsSettingsDiv,
      this.addPiechartInput.root,
      this.addFormInput.root,
    ], { style: { minWidth: '450px' }, classes: 'ui-form' });
  }
  
  public getParams() {
    return {
      table: this.tableInput.value!,
      col: this.colInput ? this.colInput.value! : null,
      templateContent: JSON.stringify(this.updatedProperties),
      models: this.tree?.items.filter((item) => !(item instanceof DG.TreeViewGroup) && item.checked).map((item) => item.text),
      addPiechart: this.addPiechartInput.value,
      addForm: this.addFormInput.value,
    };
  }

  private getTextHeight(text: string) {
    const div = ui.div();
    div.classList.add('admetica-div');
    div.textContent = text;
    document.body.appendChild(div);
    const textHeight = div.clientHeight;
    document.body.removeChild(div);
    
    return textHeight;
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