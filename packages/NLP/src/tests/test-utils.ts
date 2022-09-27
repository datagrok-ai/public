export function checkHTMLElementByInnerText(className: string, innerText: string): HTMLElement | undefined {
  const elements = document.getElementsByClassName(className);
  let check = false;
  let element;
  for (let i = 0; i < elements.length; i++) {
    if ((<HTMLElement>elements[i]).innerText == innerText)
      check = true;
    element = elements[i] as HTMLElement;
  }
  if (check == false)
    throw 'Element with innerText = "' + innerText + '" not found';
  else
    return element;
}

export function getHTMLElementByInnerText(className: string, innerText: string): HTMLElement | undefined {
  const elements = document.getElementsByClassName(className) as HTMLCollectionOf<HTMLElement>;
  for (let i = 0; i < elements.length; i++) {
    if (elements[i].innerText == innerText)
      return elements[i];
  }
}
