import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

const window = grok.shell.windows;

const topmenuToggle = ui.div([ui.iconFA('window-maximize')], 'windows-manager-toggle');
const toolboxToogle = ui.div([ui.iconFA('ballot')], 'windows-manager-toggle');
const propertiesToggle = ui.div([ui.iconFA('sliders-h')], 'windows-manager-toggle');
const helpToggle = ui.div([ui.iconFA('info')], 'windows-manager-toggle');
const vairablesToggle = ui.div([ui.iconFA('value-absolute')], 'windows-manager-toggle');
const consoleToggle = ui.div([ui.iconFA('terminal')], 'windows-manager-toggle');
const presentationToggle = ui.div([ui.iconFA('presentation')], 'windows-manager-toggle');

presentationToggle.addEventListener('click', ()=> {
  window.presentationMode ? window.presentationMode = false : window.presentationMode = true;
  setToggleState(window.presentationMode, presentationToggle);
});

toolboxToogle.addEventListener('click', ()=> {
  window.showToolbox ? window.showToolbox = false : window.showToolbox = true;
  setToggleState(window.showToolbox, toolboxToogle);
});

topmenuToggle.addEventListener('click', ()=> {
  window.simpleMode ? window.simpleMode = false : window.simpleMode = true;
  window.simpleMode ? topmenuToggle.className = 'windows-manager-toggle' : topmenuToggle.className = 'windows-manager-toggle active';
});

propertiesToggle.addEventListener('click', ()=> {
  window.showProperties ? window.showProperties = false : window.showProperties = true;
  setToggleState(window.showProperties, propertiesToggle);
});

helpToggle.addEventListener('click', ()=> {
  window.showHelp ? window.showHelp = false : window.showHelp = true;
  setToggleState(window.showHelp, helpToggle);
});

vairablesToggle.addEventListener('click', ()=> {
  window.showVariables ? window.showVariables = false : window.showVariables = true;
  setToggleState(window.showVariables, vairablesToggle);
});

consoleToggle.addEventListener('click', ()=> {
  window.showConsole ? window.showConsole = false : window.showConsole = true;
  setToggleState(window.showConsole, consoleToggle);
});

function setToggleState(v: boolean, toggle: HTMLDivElement) {
  return (v ? toggle.className = 'windows-manager-toggle active' : toggle.className = 'windows-manager-toggle');
}

function setButtonsToggleState() {
  window.simpleMode ? topmenuToggle.className = 'windows-manager-toggle' : topmenuToggle.className = 'windows-manager-toggle active';
  setToggleState(window.showToolbox, toolboxToogle);
  setToggleState(window.showProperties, propertiesToggle);
  setToggleState(window.showHelp, helpToggle);
  setToggleState(window.showVariables, vairablesToggle);
  setToggleState(window.showConsole, consoleToggle);
  setToggleState(window.presentationMode, presentationToggle);
}

export function windowsManagerPanel() {
  const root = ui.div([
    ui.tooltip.bind(topmenuToggle, () => ui.div(['Tabs ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(toolboxToogle, () => ui.div(['Toolbox ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(propertiesToggle, () => ui.div(['Context Panel ', ui.span(['F4'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(helpToggle, () => ui.div(['Context Help ', ui.span(['F1'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(vairablesToggle, () => ui.div(['Variables ', ui.span(['ALT+V'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(consoleToggle, () => ui.div(['Console ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(presentationToggle, () => ui.div(['Presentation mode ', ui.span(['F7'], {style: {color: 'var(--grey-4)'}})]), 'top'),
  ]);

  root.className = 'windows-manager-statusbar';
  document.getElementsByClassName('d4-global-status-panel')[0]?.append(root);
  setButtonsToggleState();
}

grok.events.onEvent('grok-panels-changed').subscribe((_) => setButtonsToggleState());
