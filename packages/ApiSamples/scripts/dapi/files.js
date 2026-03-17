//help-url: https://datagrok.ai/help/develop/how-to/access-data
// JS API methods for working with files

(async () => {

  let res = null;

  // Write files
  await grok.dapi.files.write('System:AppData/ApiSamples/testFile.dat', [0, 1, 2]);
  await grok.dapi.files.writeAsText('System:AppData/ApiSamples/testFile.txt', 'testString');

  // Rename
  await grok.dapi.files.writeAsText('System:AppData/ApiSamples/forRename.txt', 'testString');
  await grok.dapi.files.rename('System:AppData/ApiSamples/forRename.txt', 'renamed.txt');

  // Read files
  res = await grok.dapi.files.readAsBytes('System:AppData/ApiSamples/testFile.dat');
  console.log(`readAsBytes: ${res}`);
  res = await grok.dapi.files.readAsText('System:AppData/ApiSamples/testFile.txt');
  console.log(`readAsText: ${res}`);
  res = await grok.dapi.files.exists('System:AppData/ApiSamples/testFile.dat');
  console.log(`testFile.dat exists: ${res}`);

  // Search files
  let recursive = true;
  let searchPattern = 'test';
  res = await grok.dapi.files.list('System:AppData/ApiSamples/', recursive, searchPattern);
  console.log(`list: ${res}`);

  //Create Directory
  await grok.dapi.files.createDirectory('System:AppData/ApiSamples/geo');
  
  // Move files 
  await grok.dapi.files.move(['System:AppData/ApiSamples/testFile.txt'], 'ApiSamples/geo');
  res = await grok.dapi.files.exists('System:AppData/ApiSamples/geo/testFile.txt');
  console.log(`testFile.txt was moved to geo: ${res}`);

  //Remove Directory
  await grok.dapi.files.delete('System:AppData/ApiSamples/geo');
  
  // Delete files
  await grok.dapi.files.delete('System:AppData/ApiSamples/geo/testFile.txt');
  await grok.dapi.files.delete('System:AppData/ApiSamples/testFile.dat');
  await grok.dapi.files.delete('System:AppData/ApiSamples/renamed.txt');
})();
