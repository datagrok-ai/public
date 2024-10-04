import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

/** Max waiting time */
const MAX_TIME = 10000;

/** Timeout tick */
const TICK = 100;

/** Get the specified child */
export async function getElement(parent: HTMLElement, selectors: string) {
  let element: HTMLElement;
  let totalTime = 0;

  while (true) {
    element = parent.querySelector(selectors);

    if (element)
      break;

    if (totalTime > MAX_TIME)
      return null;

    await new Promise((resolve) => setTimeout(resolve, TICK));
    totalTime += TICK;
  }

  return element;
}

/** Get the specified child */
export async function getView(name: string) {
  let view: DG.ViewBase;
  let totalTime = 0;

  while (true) {
    view = grok.shell.v;

    if (view.name.includes(name))
      return view;

    if (totalTime > MAX_TIME)
      return null;

    await new Promise((resolve) => setTimeout(resolve, TICK));
    totalTime += TICK;
  }
}
