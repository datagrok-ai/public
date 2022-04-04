  export function updateDivInnerHTML(div: HTMLElement, content: any) {
    div.innerHTML = '';
    div.append(content);
  }