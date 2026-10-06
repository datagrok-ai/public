import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {getCurrentUserGroup} from './spotlight/group-favorites';

const window = grok.shell.windows;

/** A status-bar toggle, named and ARIA-labelled so it is addressable by what it does
 * rather than by the FontAwesome class of its icon. */
function toggle(icon: string, name: string): HTMLDivElement {
  const root = ui.div([ui.iconFA(icon)], 'windows-manager-toggle');
  root.setAttribute('name', `toggle-${name.toLowerCase().replace(/ /g, '-')}`);
  root.setAttribute('role', 'button');
  root.setAttribute('aria-label', name);
  root.setAttribute('aria-pressed', 'false');
  return root;
}

const aiToggle = toggle('user-robot', 'AI');
const topmenuToggle = toggle('window-maximize', 'Tabs');
const toolboxToogle = toggle('ballot', 'Toolbox');
const propertiesToggle = toggle('sliders-h', 'Context Panel');
const helpToggle = toggle('info', 'Context Help');
const vairablesToggle = toggle('value-absolute', 'Variables');
const consoleToggle = toggle('terminal', 'Console');
const presentationToggle = toggle('presentation', 'Presentation mode');
const inspectorToggle = toggle('tools', 'Inspector');

presentationToggle.addEventListener('click', ()=> {
  window.presentationMode ? window.presentationMode = false : window.presentationMode = true;
  setToggleState(window.presentationMode, presentationToggle);
});

toolboxToogle.addEventListener('click', ()=> {
  window.showToolbox ? window.showToolbox = false : window.showToolbox = true;
  setToggleState(window.showToolbox, toolboxToogle);
});

aiToggle.addEventListener('click', ()=> {
  window.showAI = !window.showAI;
  setToggleState(window.showAI, aiToggle);
});

topmenuToggle.addEventListener('click', ()=> {
  window.simpleMode ? window.simpleMode = false : window.simpleMode = true;
  setToggleState(!window.simpleMode, topmenuToggle);
});

propertiesToggle.addEventListener('click', ()=> {
  window.showContextPanel ? window.showContextPanel = false : window.showContextPanel = true;
  setToggleState(window.showContextPanel, propertiesToggle);
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

inspectorToggle.addEventListener('click', () => {
  if (isInspectorVisible())
    closeInspectorPane();
  else {
    document.dispatchEvent(new KeyboardEvent('keydown', {keyCode: 73, altKey: true, bubbles: true}));
    setTimeout(() => setToggleState(isInspectorVisible(), inspectorToggle), 300);
  }
});

function isInspectorVisible(): boolean {
  return Array.from(document.querySelectorAll('.tab-handle-text')).some((el) => el.textContent?.trim() === 'Inspector');
}

function closeInspectorPane(): void {
  for (const tab of Array.from(document.querySelectorAll('.tab-handle')))
    if (tab.querySelector('.tab-handle-text')?.textContent?.trim() === 'Inspector')
      (tab.querySelector('.tab-handle-close-button') as HTMLElement)?.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
}

function setToggleState(v: boolean, toggle: HTMLDivElement) {
  ui.setClass(toggle, 'active', v);
  toggle.setAttribute('aria-pressed', `${v}`);
}

function setButtonsToggleState() {
  setToggleState(window.showAI, aiToggle);
  setToggleState(!window.simpleMode, topmenuToggle);
  setToggleState(window.showToolbox, toolboxToogle);
  setToggleState(window.showContextPanel, propertiesToggle);
  setToggleState(window.showHelp, helpToggle);
  setToggleState(window.showVariables, vairablesToggle);
  setToggleState(window.showConsole, consoleToggle);
  setToggleState(window.presentationMode, presentationToggle);
  setToggleState(isInspectorVisible(), inspectorToggle);
}

export async function windowsManagerPanel() {
  const userGroup = await getCurrentUserGroup();
  const isDeveloper = userGroup?.memberships.some((g) => g.id === DG.Group.defaultGroupsIds.Developers) ?? false;

  const toggles: HTMLElement[] = [
    ui.tooltip.bind(aiToggle, () => ui.div(['AI ', ui.span(['Ctrl+I'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(topmenuToggle, () => ui.div(['Tabs ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(toolboxToogle, () => ui.div(['Toolbox ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(propertiesToggle, () => ui.div(['Context Panel ', ui.span(['F4'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(helpToggle, () => ui.div(['Context Help ', ui.span(['F1'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(vairablesToggle, () => ui.div(['Variables ', ui.span(['ALT+V'], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(consoleToggle, () => ui.div(['Console ', ui.span([''], {style: {color: 'var(--grey-4)'}})]), 'top'),
    ui.tooltip.bind(presentationToggle, () => ui.div(['Presentation mode ', ui.span(['F7'], {style: {color: 'var(--grey-4)'}})]), 'top'),
  ];
  if (isDeveloper)
    toggles.unshift(ui.tooltip.bind(inspectorToggle, () => ui.div(['Inspector ', ui.span(['ALT+I'], {style: {color: 'var(--grey-4)'}})]), 'top'));

  const root = ui.div(toggles);
  root.className = 'windows-manager-statusbar';
  document.getElementsByClassName('d4-global-status-panel')[0]?.append(root);
  setButtonsToggleState();
}

grok.events.onEvent('grok-panels-changed').subscribe((_) => setButtonsToggleState());
