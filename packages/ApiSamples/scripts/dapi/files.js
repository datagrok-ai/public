//name: Files
//description: JS API: Work with files demo
//language: javascript

// write files
grok.dapi.files.write('Demo:TestJobs:Files:DemoFiles/testFile.dat', [0,1,2]);
grok.dapi.files.writeAsText('Demo:TestJobs:Files:DemoFiles/testFile.dat/testFile.txt', 'testString');

//read files
grok.dapi.files.readAsBytes('Demo:TestJobs:Files:DemoFiles/testFile.dat');
grok.dapi.files.readAsText('Demo:TestJobs:Files:DemoFiles/testFile.txt');

grok.dapi.files.exists('Demo:TestJobs:Files:DemoFiles/testFile.dat')
    .then(isFile => console.log(isFile));

let recursive = true;
let searchPattern = "";
grok.dapi.files.list('Demo:TestJobs:Files:DemoFiles/geo', recursive, searchPattern);
grok.dapi.files.move(['Demo:TestJobs:Files:DemoFiles/testFile.txt'],'/geo');

grok.dapi.files.delete('Demo:TestJobs:Files:DemoFiles/testFile.dat/testFile.txt');


