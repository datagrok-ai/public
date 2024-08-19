//tags: FileInfo
//help-url: https://datagrok.ai/help/develop/how-to/access-data
// JS API methods for working with files

(async () => {

  let res = null;

  // Write files
  await grok.dapi.files.write('Samples:Files/testFile.dat', [0, 1, 2]);
  await grok.dapi.files.writeAsText('Samples:Files/testFile.txt', 'testString');

  // Rename
  await grok.dapi.files.writeAsText('Samples:Files/forRename.txt', 'testString');
  await grok.dapi.files.rename('Samples:Files/forRename.txt', 'renamed.txt');

  // Read files
  res = await grok.dapi.files.readAsBytes('Samples:Files/testFile.dat');
  console.log(`readAsBytes: ${res}`);
  res = await grok.dapi.files.readAsText('Samples:Files/testFile.txt');
  console.log(`readAsText: ${res}`);
  res = await grok.dapi.files.exists('Samples:Files/testFile.dat');
  console.log(`testFile.dat exists: ${res}`);

  // Search files
  let recursive = true;
  let searchPattern = 'world';
  res = await grok.dapi.files.list('Samples:Files/geo', recursive, searchPattern);
  console.log(`list: ${res}`);

  // Move files
  await grok.dapi.files.move(['Samples:Files/testFile.txt'], 'geo');
  res = await grok.dapi.files.exists('Samples:Files/geo/testFile.txt');
  console.log(`testFile.txt was moved to geo: ${res}`);

  // Delete files
  await grok.dapi.files.delete('Samples:Files/geo/testFile.txt');
  await grok.dapi.files.delete('Samples:Files/testFile.dat');
})();
