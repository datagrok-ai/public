import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import '../../css/chem.css';

export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export function addCopyIcon(element: Object | HTMLElement, paneName: string) {
  const copyIcon = ui.icons.copy(() => {}, 'Copy');
  copyIcon.onclick = (e) => {
    let text: string;
    if (element instanceof HTMLElement) {
      text = element.innerText;
    } else {
      let tableString = '';
      for (const [key, value] of Object.entries(element))
        tableString += `${key}\t${(value as any).innerText}\n`;
      text = tableString;
    }
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