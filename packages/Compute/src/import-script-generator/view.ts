import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { ImportGeneratorStore } from './store';

export class ImportScriptGeneratorApp extends DG.ViewBase {
  store = new ImportGeneratorStore();
  scriptNameInput = ui.stringInput('Script name', '', (v: string) => this.store.scriptName = v);
  descriptionInput = ui.stringInput('Description', '', (v: string) => this.store.description = v);

  constructor() {
    super();

    const runButton = ui.bigButton('Generate script', async () => await this.generate());
    const saveButton = ui.bigButton('Save', async () => await this.save());
    const openButton = ui.bigButton('Open', async () => this.open());
    openButton.style.display = 'none';

    this.store.lastSavedScriptId.subscribe((newId) => {
      if (newId !== '') {
        openButton.style.removeProperty('display');
      }
    })

    const dynamicDFInput = (counter: number) => {
      const deleteBtn = ui.button('Remove', () => {});
      deleteBtn.style.width = deleteBtn.style.height;
      deleteBtn.style.alignSelf = 'center';

      this.store.dfSheetsToSeek.push('');

      const containter =  ui.divH([
        ui.stringInput(`Sheet name #${counter}`, '', (v: string) => this.store.dfSheetsToSeek[counter] = v).root,
        deleteBtn,
      ], {style: {'justify-content': 'space-between'}});
      deleteBtn.addEventListener('click', () => {
        this.store.dfSheetsToSeek.splice(counter, 1);
        containter.remove();
      });

      return containter;
    }

    const dynamicScalarInput = (counter: number) => {
      const deleteBtn = ui.button('Remove', () => {});
      deleteBtn.style.width = deleteBtn.style.height;
      deleteBtn.style.alignSelf = 'center';

      this.store.scalarsToSeek.push('');

      const containter =  ui.divH([
        ui.stringInput(`Scalar name #${counter}`, '', (v: string) => this.store.scalarsToSeek[counter] = v).root,
        deleteBtn,
      ], {style: {'justify-content': 'space-between'}});
      deleteBtn.addEventListener('click', () => {
        this.store.scalarsToSeek.splice(counter, 1);
        containter.remove();
      });

      return containter;
    }

    const dfInputsDiv = ui.divV([]);
    const dfInputAddBtn = ui.button('Add DF to search', () => {
      dfInputsDiv.append(dynamicDFInput(this.store.dfSheetsToSeek.length))
    });

    const scalarsInputsDiv = ui.divV([]);
    const scalarsInputAddBtn = ui.button('Add scalar to search', () => {
      scalarsInputsDiv.append(dynamicScalarInput(this.store.scalarsToSeek.length))
    });

    const formDiv = ui.divV([
      ui.h3('Script properties'),
      this.scriptNameInput.root,
      this.descriptionInput.root,
      ui.h3('Dataframes to search'),
      dfInputsDiv,
      dfInputAddBtn,
      ui.h3('Scalars to search'),
      scalarsInputsDiv,
      scalarsInputAddBtn,
      ui.divH([runButton], {style: {'justify-content': 'flex-end'}})
    ], 'ui-form');

    let scriptDiv = ui.divText(this.store.scriptText.value);
    let infoDiv = ui.info(this.store.infoText.value);

    const outputDiv = ui.panel([
      ui.h3('Script text'),
      scriptDiv,
      infoDiv,
      ui.divH([saveButton, openButton], {style: {'justify-content': 'flex-end'}})
    ], {style: {display: 'none'}})
    this.store.scriptText.subscribe((newScriptText) => {
      const newScriptDiv = ui.divText(newScriptText);
      scriptDiv.replaceWith(newScriptDiv);
      scriptDiv = newScriptDiv;

      if (newScriptText !== '') outputDiv.style.removeProperty('display');
    })
    this.store.infoText.subscribe((newInfoText) => {
      const newInfoDiv = (newInfoText === '') ? ui.div([]) : ui.info(newInfoText);
      infoDiv.replaceWith(newInfoDiv);
      infoDiv = newInfoDiv;
    })

    this.root.appendChild(ui.splitH([formDiv, outputDiv], {style: {height: '100%', width: '100%'}}, true));

    ui.tools.handleResize(formDiv.parentElement!, (width) => {
      if (width < 350)
        formDiv.classList.add('ui-form-condensed');
      else
        formDiv.classList.remove('ui-form-condensed');
    });
    formDiv.parentElement!.style.maxWidth = '340px';
  } 
  
  async generate() {
    const dfOutputsLines = this.store.dfSheetsToSeek.map((dfName) => `//output: dataframe ${dfName.trim()} { viewer: Grid(); }`).join("\n");
    const scalarsLines = this.store.scalarsToSeek.map((scalarName) => `//output: string ${scalarName.trim()}`).join("\n");

    const fileParsingCode = 
`
const wb = new ExcelJS.Workbook();
await wb.xlsx.load(await uploadedFile.arrayBuffer());
const dfSheetsToSeek = ["${this.store.dfSheetsToSeek.join(`","`)}"];`

    const dfParsingCode = (dfName: string) => 
`
const ${dfName}_ws = wb.getWorksheet("${dfName}");
let ${dfName}_data = await ${dfName}_ws.workbook.csv.writeBuffer();
${dfName} = grok.data.parseCsv(${dfName}_data.toString());`

    const getAllScalars = () => 
`
const scalarSheets = wb.worksheets.filter(ws => !dfSheetsToSeek.includes(ws.name));
const scalarNames = scalarSheets.flatMap((sheet) => sheet.getColumn(1).values);
const scalarValues = scalarSheets.flatMap((sheet) => sheet.getColumn(2).values);`

    const scalarParsingCode = (searchName: string) => 
`
const ${searchName}_index = scalarNames.findIndex((scalarName) => scalarName === "${searchName}");
${searchName} = (${searchName}_index >= 0) ? scalarValues[${searchName}_index] : '';`
    
    let scriptText = `//name: ${this.store.scriptName}` + "\n";
    scriptText += `//description: ${this.store.description}` + "\n";
    scriptText += `//language: javascript` + "\n";
    scriptText += `//input: file uploadedFile` + "\n";
    scriptText += scalarsLines + "\n";
    scriptText += dfOutputsLines + "\n";
    scriptText += `//editor: Compute:FunctionViewEditor` + "\n";
    scriptText += fileParsingCode + "\n";

    this.store.dfSheetsToSeek.forEach((dfName) => {
      scriptText += dfParsingCode(dfName.trim()) + "\n";
    });

    scriptText += getAllScalars();

    this.store.scalarsToSeek.forEach((scalarName) => {
      scriptText += scalarParsingCode(scalarName) + "\n";
    })

    this.store.scriptText.next(scriptText);

    const duplicates = await grok.dapi.scripts.filter(`name="${this.store.scriptName}"`).list();
    if (duplicates.length > 0){
      this.store.infoText.next(`Script named "${this.store.scriptName}" already exists`);
    }else {
      this.store.infoText.next(``);
    }
  }

  async save() {
    const resultScript = DG.Script.create(this.store.scriptText.value);
    const savedScript = await grok.dapi.scripts.save(resultScript);
    this.store.lastSavedScriptId.next(savedScript.id);
    grok.shell.info(`Script "${this.store.scriptName}" has been saved successfully`)
  }

  open() {
    window.open(`${window.location.origin}/script/${this.store.lastSavedScriptId.value}`)
  }
}