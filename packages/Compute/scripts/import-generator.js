//name: Import script generator
//language: javascript
//tags: higher-order function
//input: string dfSheetsToSeek
//input: string scalarSheetsToSeek
//input: string scalarsToSeek 
//input: string scriptName 
//input: string description 

const dfOutputsLines = dfSheetsToSeek.split(',').map((dfName) => `//output: dataframe ${dfName.trim()} { viewer: Grid(); }`).join("\n");
const scalarsLines = scalarsToSeek.split(',').map((scalarName) => `//output: string ${scalarName.trim()}`).join("\n");

let scriptText = `//name: ${scriptName}` + "\n";
scriptText += `//description: ${description}` + "\n";
scriptText += `//language: javascript` + "\n";
scriptText += `//input: string file {type: file}` + "\n";
scriptText += scalarsLines + "\n";
scriptText += dfOutputsLines + "\n";
scriptText += `//editor: Compute:FunctionViewEditor` + "\n";
scriptText += 

const duplicates = await grok.dapi.scripts.filter(`name="${scriptName}"`).list();

const dialog = ui.dialog('Confirm script text');
dialog.add(ui.h3('Script text'));
dialog.add(ui.divText(scriptText));

if (duplicates.length > 0) {
    dialog.add(ui.info(`Script named ${scriptName} already exists`));
}

dialog.onOK(async () => {
    const resultScript = DG.Script.create(scriptText);
    await grok.dapi.scripts.save(resultScript);
    grok.shell.info(`Script "${scriptName}" has been saved successfully`)
})
dialog.show({center: true});  