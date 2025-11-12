import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askAiTableView} from './viewers-ai';
import {IDartApi} from "datagrok-api/src/api/grok_api.g";

export const aiPanel: HTMLElement = ui.divV([]);
const input: HTMLTextAreaElement = ui.element('textarea');
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export function initAiPanel() {
  input.placeholder = 'Ask a question, or tell what to do...';
  input.style.overflow = 'none';
  input.style.boxShadow = 'none';
  input.style.borderColor = '--var(grey-2)'
  input.style.width = '100%';
  const spinner = ui.icons.spinner();

  const prompt = async (question: string) => {
    aiPanel.append(ui.divText(question, {classes: 'ai-panel-question'}));

    if (api.grok_Shell_AI_Prompt(question)) {
      return;
    }

    const currentView = grok.shell.v;
    if (currentView && currentView.type === DG.VIEW_TYPE.TABLE_VIEW && question) {
      const tableView = currentView as DG.TableView;
      aiPanel.append(spinner);
      try {
        await askAiTableView(tableView, input.value);
      }
      finally {
        spinner.remove();
      }
      input.value = '';
    }
  }

  // Add Enter key event listener
  input.addEventListener('keydown', async (event) => {
    if (event.key === 'Enter' && !event.shiftKey) {
      event.preventDefault();
      prompt(input.value.trim());
    }
  });

  const mikeIcon = ui.iconFA('microphone', () => {
    const sr = new SpeechRecognition();
    sr.lang = 'en-US';
    sr.onresult = (x) => {
      grok.shell.info(x);
      const transcript = x.results[0][0].transcript;
      prompt(transcript);
    };
    sr.start();
  });

  aiPanel.append(ui.divH([mikeIcon]));
  aiPanel.append();
  aiPanel.append(input);
}

export function setAiPanelVisibility(visible: boolean) {
  if (visible) {
    const contextRoot = grok.shell.windows.context.root;
    const contextNode = contextRoot ? grok.shell.dockManager.findNode(contextRoot) : null;
    //@ts-ignore
    grok.shell.dockManager.dock(aiPanel, contextNode ? 'down' : 'right', contextNode, 'AI', 0.25);
    input.tabIndex = 0;
    input.focus();
  }
  else {
    grok.shell.dockManager.close(aiPanel);
  }
}

export function getAiPanelVisibility(): boolean {
  return document.body.contains(aiPanel);
}