export function updateDivInnerHTML(div: HTMLElement, content: string | Node): void {
  div.innerHTML = '';
  div.append(content);
}
