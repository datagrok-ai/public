//name: Files
//description: JS API: Work with files demo
//language: javascript

(async () => {

let res = null;

// write files
await grok.dapi.files.write('Demo:TestJobs:Files:DemoFiles/testFile.dat', [0, 1, 2]);
await grok.dapi.files.writeAsText('Demo:TestJobs:Files:DemoFiles/testFile.txt', 'testString');

//rename
await grok.dapi.files.writeAsText('Demo:TestJobs:Files:DemoFiles/forRename.txt', 'testString');
await grok.dapi.files.rename('Demo:TestJobs:Files:DemoFiles/forRename.txt', 'renamed.txt');

//read files

res = grok.dapi.files.readAsBytes('Demo:TestJobs:Files:DemoFiles/testFile.dat');
console.log(`readAsBytes: ${res}`);
res = grok.dapi.files.readAsText('Demo:TestJobs:Files:DemoFiles/testFile.txt');
console.log(`readAsText: ${res}`);
res = grok.dapi.files.exists('Demo:TestJobs:Files:DemoFiles/testFile.dat');
console.log(`exists: ${res}`);

let recursive = true;
let searchPattern = "";
res = await grok.dapi.files.list('Demo:TestJobs:Files:DemoFiles/geo', recursive, searchPattern);
console.log(`list: ${res}`);

await grok.dapi.files.move(['Demo:TestJobs:Files:DemoFiles/testFile.txt'],'/geo');

await grok.dapi.files.delete('Demo:TestJobs:Files:DemoFiles/geo/testFile.txt');
await grok.dapi.files.delete('Demo:TestJobs:Files:DemoFiles/testFile.dat');

})();