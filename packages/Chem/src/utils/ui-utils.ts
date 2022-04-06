import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

export function updateDivInnerHTML(div: HTMLElement, content: any) {
  div.innerHTML = '';
  div.append(content);
}


export function updateMetricsLink(distance: string, fingerprint: string, metricsDiv: HTMLDivElement, object: any, options: any) {
  const metricsButton = ui.button(`${distance}/${fingerprint}`, () => {
    if (!grok.shell.windows.showProperties)
      grok.shell.windows.showProperties = true;
    grok.shell.o = object;
  });
  //@ts-ignore
  Object.keys(options).forEach((it) => metricsButton.style[it] = options[it]);
  updateDivInnerHTML(metricsDiv, metricsButton);
}
