//Two column layout

let windows = grok.shell.windows;
windows.showToolbox = false; //*hide tool box
windows.showContextPanel = false; //*hide context panel
windows.showHelp = false; //*hide context help

let col2 = ui.block50([]);
let col3 = ui.block([]);

let scripts = ui.divH([], {style:{flexWrap:'wrap'}});

let scriptViewer = ui.textInput('','');
scriptViewer.input.style.height = '300px';
scriptViewer.readOnly = true;
scriptViewer.input.style.resize = 'none';
console.log(scriptViewer);

let scriptDetails = ui.div();


(async()=>{
    await grok.dapi.scripts.list().then((script)=>
      script.map((p)=> {
      let scriptData = ui.renderCard(p);
      scriptData.addEventListener('click', function (event) {
       scriptViewer.value = p.script;
       grok.shell.info(p.name);
       scriptDetails.innerHTML = '';
       scriptDetails.append(ui.div([
         ui.tableFromMap({
           'Name': p.name,
           'Description': p.description
         }),
         ui.p(),
         ui.bigButton('Run', ()=>{p.apply()}),
         ui.link('Open script',p.path,'')
       ], {style:{marginTop:'15px'}})) 
      })
      scriptData.classList += ' ui-div ui-block ui-block-50';
      scriptData.style.margin = '0px';
      scriptData.style.display = 'flex';
      scriptData.style.alignItems = 'strech';
      scripts.append(scriptData);
      
      
    }));
    })();


col2.append(ui.panel([
  ui.h1('Scripts'),
  scripts,
]));

col3.append(ui.panel([
  ui.h1('Script viewer'),
  ui.block(ui.box(scriptViewer.root)),
  scriptDetails
]));

let view = grok.shell.newView('Layout 2', [col2,col3]);
view.box = true;