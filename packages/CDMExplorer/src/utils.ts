export function sortObject(o) {
    return Object.keys(o).sort().reduce((r, k) => (r[k] = o[k], r), {});
}

export function updateDivInnerHTML(div: HTMLElement, content: any) {
    div.innerHTML = '';
    div.append(content);
  }
