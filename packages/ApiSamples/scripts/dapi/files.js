//name: Files
//description: JS API: Work with files demo
//language: javascript

// write files
grok.dapi.files.write('Demo:TestJobs:Files:DemoFiles/testFile.dat', [0,1,2]);
grok.dapi.files.writeAsText('Demo:TestJobs:Files:DemoFiles/testFile.txt', 'testString');

//rename
grok.dapi.files.writeAsText('Demo:TestJobs:Files:DemoFiles/forRename.txt', 'testString')
    .then(grok.dapi.files.rename('Demo:TestJobs:Files:DemoFiles/forRename.txt', 'renamed.txt'));

//read files
grok.dapi.files.readAsBytes('Demo:TestJobs:Files:DemoFiles/testFile.dat').then(res => console.log(res));
grok.dapi.files.readAsText('Demo:TestJobs:Files:DemoFiles/testFile.txt').then(res => console.log(res));

grok.dapi.files.exists('Demo:TestJobs:Files:DemoFiles/testFile.dat').then(isFile => console.log(isFile));

let recursive = true;
let searchPattern = "";
grok.dapi.files.list('Demo:TestJobs:Files:DemoFiles/geo', recursive, searchPattern).then(res => console.log(res));;

grok.dapi.files.move(['Demo:TestJobs:Files:DemoFiles/testFile.txt'],'/geo');

grok.dapi.files.delete('Demo:TestJobs:Files:DemoFiles/geo/testFile.txt');
grok.dapi.files.delete('Demo:TestJobs:Files:DemoFiles/testFile.dat');