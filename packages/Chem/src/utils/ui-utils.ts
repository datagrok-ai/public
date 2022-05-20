export function updateDivInnerHTML(div: HTMLElement, content: any): void {
  div.innerHTML = '';
  div.append(content);
}
