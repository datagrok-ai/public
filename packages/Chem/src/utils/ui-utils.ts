import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import '../../css/chem.css';

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export function addCopyIcon(text: string | HTMLElement, paneName: string) {
  const copyIcon = ui.icons.copy(() => {}, 'Copy');
  copyIcon.onclick = (e) => {
    if (typeof(text) != 'string') 
      navigator.clipboard.writeText(text.innerText);
    else 
      navigator.clipboard.writeText(text);
    grok.shell.info('Copied to clipboard');
    e.stopImmediatePropagation();
  } 
  copyIcon.classList.add('copy-icon');
  const accPanes = document.getElementsByClassName('d4-accordion-pane-header');
  for (let i = 0; i < accPanes.length; ++i) {
    if (accPanes[i].innerHTML === paneName) {
      const pane = accPanes[i];
      pane.append(copyIcon);
      pane.parentElement?.addEventListener('mouseenter', () => {copyIcon.style.visibility = 'visible'});
      pane.parentElement?.addEventListener('mouseleave', () => {copyIcon.style.visibility = 'hidden'});
    }
  }
}