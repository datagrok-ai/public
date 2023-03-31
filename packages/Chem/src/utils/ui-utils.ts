export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}

export const copyIconStyles = {
  position: 'absolute', 
  left: '170px', 
  top: '-20px'
}