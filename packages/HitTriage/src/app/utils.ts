import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subscription} from 'rxjs';

export const toFormatedDateString = (d: Date): string => {
  return `${d.getFullYear()}/${d.getMonth() + 1}/${d.getDate()}`;
};

export function modifyUrl(key: string, value: string) {
  const title = document.title;
  const url = window.location.href.split('?')[0] + '?' + key + '=' + value;
  if (history.replaceState) {
    const obj = {
      Title: title,
      Url: url,
    };
    history.replaceState(obj, obj.Title, obj.Url);
  }
}

export function hideComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'none';
  });
}

export function showComponents(styles: CSSStyleDeclaration[]) {
  styles.forEach((style) => {
    style.display = 'flex';
  });
}

export function addBreadCrumbsToRibbons(view: DG.ViewBase, firstPath: string, secondPath: string, handler: Function):
  {breadcrumbs: DG.Breadcrumbs, sub: Subscription} {
  const breadcrumbs = ui.breadcrumbs([firstPath, secondPath]);
  //breadcrumbs.root.id = 'hit-triage-breadcrumbs';
  const sub = breadcrumbs.onPathClick.subscribe((path) => {
    if (path.length === 1) {
      sub.unsubscribe();
      popRibbonPannels(view);
      handler();
    }
  });
  const viewNameRoot = grok.shell.v.ribbonMenu.root.parentElement?.getElementsByClassName('d4-ribbon-name')[0];
  if (viewNameRoot) {
    viewNameRoot.textContent = '';
    viewNameRoot.appendChild(breadcrumbs.root);
  }
  /*
  let ribbons = view.getRibbonPanels();
  if (!ribbons?.length)
    ribbons = [];
  ribbons.push([breadcrumbs.root]);

  if (document.getElementById('hit-triage-breadcrumbs') === null) {
    view.setRibbonPanels(ribbons);
    breadcrumbs.root.parentElement?.classList.add('d4-ribbon-name');
    breadcrumbs.root.parentElement?.classList.remove('d4-ribbon-item');
  }*/
  return {breadcrumbs, sub};
}

export function popRibbonPannels(view: DG.ViewBase) {
  const ribbons = view.getRibbonPanels();
  if (!ribbons?.length)
    return;
  ribbons.pop();
  view.setRibbonPanels(ribbons);
};

export function checkRibbonsHaveSubmit(ribbons: HTMLElement[][]): boolean {
  return ribbons.reduce((prev, cur) => {
    return prev || cur.reduce((p, c) =>
      p || c.classList.contains('hit-design-submit-button') ||
        Array.from(c.children).some((child) =>child.classList.contains('hit-design-submit-button')), false);
  }, false);
}
