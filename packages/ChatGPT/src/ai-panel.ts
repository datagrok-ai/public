import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {askAiTableView} from './viewers-ai';

export const aiPanel: HTMLElement = ui.divV([]);
const input: HTMLTextAreaElement = ui.element('textarea');

export function initAiPanel() {
  input.placeholder = 'Ask a question, or tell what to do...';
  input.style.overflow = 'none';
  input.style.boxShadow = 'none';
  input.style.borderColor = '--var(grey-2)'
  input.style.width = '100%';
  
  // Add Enter key event listener
  input.addEventListener('keydown', (event) => {
    if (event.key === 'Enter' && !event.shiftKey) {
      event.preventDefault();
      const currentView = grok.shell.v;
      if (currentView && currentView.type === DG.VIEW_TYPE.TABLE_VIEW && input.value.trim()) {
        const tableView = currentView as DG.TableView;
        askAiTableView(tableView, input.value);
        input.value = '';
      }
    }
  });
  
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