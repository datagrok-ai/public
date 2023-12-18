//tags: Script, Files
//help-url: https://datagrok.ai/help/compute/scripting
// An example of using files in a script

async function logFileContent(fileInfo) {
    let data = await fileInfo.readAsString();
    grok.shell.info(data);
} 

var fileInfo = DG.FileInfo.fromString('test of a script using dfiles');
logFileContent(fileInfo);