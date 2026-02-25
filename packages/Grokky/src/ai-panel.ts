import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askAiTableView} from './viewers-ai';
import {popupMenu} from 'datagrok-api/ui';

export const aiPanel: HTMLElement = ui.divV([]);
const input: HTMLTextAreaElement = ui.element('textarea');
const api: any = (typeof window !== 'undefined' ? window : global.window) as any;


export function initAiPanel() {
  input.placeholder = 'Ask a question, or tell what to do...';
  input.style.overflow = 'none';
  input.style.boxShadow = 'none';
  input.style.outline = 'none';
  input.style.borderColor = '--var(grey-2)';
  input.style.width = '100%';
  const spinner = ui.icons.spinner();
  const messagesDiv = ui.div();

  const prompt = async (question: string) => {
    messagesDiv.append(ui.divText(question, {classes: 'ai-panel-question'}));

    //@ts-ignore
    if (await api.grok_Shell_AI_Prompt(question))
      return;


    const currentView = grok.shell.v;
    if (currentView && currentView.type === DG.VIEW_TYPE.TABLE_VIEW && question) {
      const tableView = currentView as DG.TableView;
      messagesDiv.append(spinner);
      try {
        await askAiTableView(tableView, input.value);
      } finally {
        spinner.remove();
      }
      input.value = '';
    }
  };

  // Add Enter key event listener
  input.addEventListener('keydown', async (event) => {
    if (event.code === 'Enter' && !event.shiftKey) {
      console.log('prompt');
      event.preventDefault();
      prompt(input.value.trim());
      input.value = '';
    }
  });

  let listening: boolean = false;
  const recognitions: Set<SpeechRecognition> = new Set();

  function refreshUI() {
    mikeIcon.style.color = listening ? 'var(--blue-1)' : 'initial';
    listeningLabel.style.display = listening ? 'initial' : 'none';
  }

  function listen() {
    const sr = new SpeechRecognition();
    recognitions.add(sr);
    sr.lang = 'en-US';
    sr.onresult = (x) => {
      const transcript = x.results[0][0].transcript;

      if (transcript === 'stop') {
        listening = false;
        refreshUI();
      } else
        prompt(transcript);
    };
    sr.onend = () => {
      recognitions.delete(sr);
      if (listening)
        listen();
    };
    sr.start();
  }

  const mikeIcon = ui.iconFA('microphone', () => {}, 'Toggle microphone input');
  const listeningLabel = ui.divText('Listening...', {style: {
    display: 'none',
    verticalAlign: 'middle'
  }});

  const prompts = {
    'Scatter plot': [
      'set x to age',
      'y by height',
      'set title to Foo',
      'invert y axis',
      'set marker default size to 30',
      'increase marker default size',
      'show regression line',
      'hide color selector',
      'decrease marker opacity',
      'color by race',
      'set marker color to age',
      'size-code by age',
      'marker by sex',
      'zoom in',
      'select all rows',
    ],
    'Grid': [
      'next page', 'next page', 'scroll down',
      'page up', 'page down',
      'move right', 'next column',
      'scroll to the top', 'scroll to the bottom',
      'select column age', 'select column weight',
      'color-code height',
      'sort by age',
    ],
    'Table view': [
      'activate spreadsheet',
      'focus on scatter plot',
      'spreadsheet',
      'grid',
      'table',
    ],
    'Shell': [
      'open plugins view',
      'show users',
      'show plugins'
    ]
  };

  const convosIcon = ui.icons.settings(() => {
    const menu = ui.popupMenu();
    for (const name of Object.keys(prompts))
      menu.group(name).items(prompts[name as keyof typeof prompts], (s) => prompt(s));

    menu.show();
  });

  const toggle = () => {
    listening = !listening;
    refreshUI();

    if (listening)
      listen();
    else {
      for (const sr of recognitions.values())
        sr.stop();
    }
  };

  mikeIcon.onclick = () => toggle();
  listeningLabel.onclick = () => toggle();

  const clearIcon = ui.iconFA('trash', () => ui.empty(messagesDiv), 'Clear console');

  aiPanel.append(ui.divH([clearIcon, mikeIcon, convosIcon, listeningLabel], 'd4-ribbon-item'));
  aiPanel.append();
  aiPanel.append(input);
  aiPanel.append(messagesDiv);
}

export function setAiPanelVisibility(visible: boolean) {
  if (visible) {
    const contextRoot = grok.shell.windows.context.root;
    const contextNode = contextRoot ? grok.shell.dockManager.findNode(contextRoot) : null;
    //@ts-ignore
    grok.shell.dockManager.dock(aiPanel, contextNode ? 'down' : 'right', contextNode, '', 0.4);
    input.tabIndex = 0;
    input.focus();
  } else
    grok.shell.dockManager.close(aiPanel);
}

export function getAiPanelVisibility(): boolean {
  return document.body.contains(aiPanel);
}
