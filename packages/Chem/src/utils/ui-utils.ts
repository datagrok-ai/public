import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import '../../css/chem.css';

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export function addCopyIcon(text: string, paneName: string) {
  const copyIcon = ui.icons.copy(() => {
    navigator.clipboard.writeText(text);
    grok.shell.info('Copied to clipboard');
  }, 'Copy');
  
  copyIcon.classList.add('copy-icon');

  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === paneName) {
      const pane = accPanes[i];
      pane.append(copyIcon);
      pane.addEventListener('mouseenter', () => {copyIcon.style.visibility = 'visible'});
      pane.addEventListener('mouseleave', () => {copyIcon.style.visibility = 'hidden'});
    }
  }
}