//help-url: https://datagrok.ai/help/compute/scripting
// An example of using files in a script

async function logFileContentAsString(fileInfo) {
  let data = await fileInfo.readAsString();
  grok.shell.info(data);
}

async function logFileContentAsBytes(fileInfo) {
  let data = await fileInfo.readAsBytes();
  grok.shell.info(data);
}

// these values are the same for utf-8
let line = 'test how script uses files';
let fileName = 'text.txt';
let bytes =	new Uint8Array([0x74, 0x65, 0x73, 0x74, 0x20, 0x68, 0x6F, 0x77, 0x20, 0x73, 0x63, 0x72, 0x69, 0x70, 0x74, 0x20, 0x75, 0x73, 0x65, 0x73, 0x20, 0x66, 0x69, 0x6C, 0x65, 0x73]);
let fileInfo = DG.FileInfo.fromString(fileName, line);
logFileContentAsString(fileInfo);
logFileContentAsBytes(fileInfo);
fileInfo = DG.FileInfo.fromBytes(fileName, bytes);
logFileContentAsString(fileInfo);
logFileContentAsBytes(fileInfo);
