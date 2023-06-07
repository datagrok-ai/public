import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

let propertiesToggle = ui.div([
    ui.iconFA('chevron-left')
], 'windows-manager-toggle');

let helpToggle = ui.div([
    ui.iconFA('info')
], 'windows-manager-toggle');

let vairablesToggle = ui.div([
    ui.iconFA('info')
], 'windows-manager-toggle');

let consoleToggle = ui.div([
    ui.iconFA('info')
], 'windows-manager-toggle');

let windows = [
    { 
     'state': grok.shell.windows.showProperties,
     'root': propertiesToggle 
    },
    {
        'state': grok.shell.windows.showHelp,
        'root': helpToggle
    },
    {
        'state': grok.shell.windows.showVariables,
        'root': vairablesToggle
    },
    {
        'state': grok.shell.windows.showConsole,
        'root': consoleToggle
    }
];

export function windowsSidebar(){

    propertiesToggle.addEventListener('click', ()=> {
        if (grok.shell.windows.showProperties) {
            propertiesToggle.classList.remove('active');
            grok.shell.windows.showProperties = false;
        } else {
            propertiesToggle.classList.add('active');
            grok.shell.windows.showProperties = true;
        }
    });

    helpToggle.addEventListener('click', ()=> {
        if (grok.shell.windows.showHelp) {
            helpToggle.classList.remove('active');
            grok.shell.windows.showHelp = false;
        } else {
            helpToggle.classList.add('active');
            grok.shell.windows.showHelp = true;
        }
    });

    vairablesToggle.addEventListener('click', ()=> {
        if (grok.shell.windows.showVariables) {
            vairablesToggle.classList.remove('active');
            grok.shell.windows.showVariables = false;
        } else {
            vairablesToggle.classList.add('active');
            grok.shell.windows.showVariables = true;
        }
    });

    consoleToggle.addEventListener('click', ()=> {
        if (grok.shell.windows.showConsole) {
            consoleToggle.classList.remove('active');
            grok.shell.windows.showConsole = false;
        } else {
            consoleToggle.classList.add('active');
            grok.shell.windows.showConsole = true;
        }
    })

    for (let i=0; i<windows.length; i++){
        if(windows[i].state){
            windows[i].root.classList.add('active')
        }else{
            windows[i].root.classList.remove('active')
        }
    }   

    let root = ui.div([propertiesToggle, helpToggle, vairablesToggle, consoleToggle]);
    root.className = 'windows-manager-sidebar';

    document.getElementsByClassName('layout-toolbox-workarea-wrapper')[0].append(root);
}

export function windowsToolbar(){

}