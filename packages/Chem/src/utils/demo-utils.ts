import { delay } from "@datagrok-libraries/utils/src/test";

export const demoScaffold = `
Actelion Java MolfileCreator 1.0

  7  7  0  0  0  0  0  0  0  0999 V2000
   11.3750  -10.5625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.3750  -12.0625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.6740  -12.8125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.9731  -12.0625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.9731  -10.5625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.6740   -9.8125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.6740   -8.3125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  1  1  0  0  0  0
  6  7  1  0  0  0  0
M  END
`


export async function scrollTable(el: HTMLElement, delta: number, cycles: number, secDelay: number) {
    for (let i = 0; i < cycles; i++) {
        el.dispatchEvent(new WheelEvent('wheel', { deltaY: delta }));
        await delay(secDelay);
    }
}

export async function openSketcher(parent: HTMLElement) {
    const sketcherFilter = parent.getElementsByClassName('sketch-link')[0];
    sketcherFilter!.classList.add('chem-demo-sketcher-link');
    await delay(1000);
    (sketcherFilter as HTMLElement).click();
    sketcherFilter!.classList.remove('chem-demo-sketcher-link');
    await delay(1000);
    return document.getElementsByClassName('d4-dialog')[0];
}


export function getAccordionPane(name: string, parent: Element) {
    const paneHeader = Array.from(parent!.getElementsByClassName('d4-accordion-pane-header'))
        .find(el => el.textContent === name) as HTMLElement;
    if (!paneHeader!.classList.contains('expanded'))
        paneHeader!.click();
    return paneHeader.parentElement?.children[1];
}

export function closeAllAccordionPanes(parent: Element) {
    const paneHeaders = Array.from(parent!.getElementsByClassName('d4-accordion-pane-header'));
    paneHeaders.forEach((panelHeader) => {
        if (panelHeader!.classList.contains('expanded'))
            (panelHeader as HTMLElement)!.click();
    })
}
