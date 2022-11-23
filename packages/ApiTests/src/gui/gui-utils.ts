import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


export function isExceptionElement(action: string): void {
  const exceptionElement = document.getElementsByClassName('grok-global-exception');
  if (exceptionElement.length != 0)
    throw new Error('exception was thrown after action: ' + action);
}

export function isColumnPresent(columns: DG.ColumnList, columnName: string): void {
  let check = false;
  if (columns.byName(columnName) == null)
    check = true;
  if (check == true)
    throw new Error('column ' + columnName + ' not found');
}

export function isViewerPresent(viewers: DG.Viewer[], viewerName: string): void {
  let check = false;
  for (let i:number = 0; i < viewers.length; i++) {
    if (viewers[i].type == viewerName) {
      check = true;
      break;
    }
  }
  if (check == false)
    throw new Error(viewerName + ' not found');
}

export function checkViewer(viewers: DG.Viewer[], viewerName: string): boolean {
  return viewers.filter((v) => v.type == viewerName).length > 0;
}

export function isErrorBallon(ballonText: string): void {
  const exceptionElement = document.getElementsByClassName('d4-balloon-content')[0] as HTMLElement;
  if (exceptionElement == undefined)
    throw new Error('Error ballon didn\'t appear where it should have.');
  if (exceptionElement.textContent != ballonText)
    throw new Error('Wrong error message on balloon');
  exceptionElement.click();
}

export function isDialogPresent(dialogTitle: string): boolean {
  let check = false;
  for (let i=0; i < DG.Dialog.getOpenDialogs().length; i++) {
    if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
      check = true;
      break;
    }
  }
  if (check == false)
    throw new Error('Dialog ' + dialogTitle + ' not found');
  return check;  
}

export function checkDialog(dialogTitle: string): boolean {
  let check = false;
  for (let i=0; i < DG.Dialog.getOpenDialogs().length; i++) {
    if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
      check = true;
      break;
    }
  }
  return check;  
}

export function returnDialog(dialogTitle: string): DG.Dialog | undefined {
  let dialog: DG.Dialog | undefined;
  for (let i = 0; i < DG.Dialog.getOpenDialogs().length; i++) {
    if (DG.Dialog.getOpenDialogs()[i].title == dialogTitle) {
      dialog = DG.Dialog.getOpenDialogs()[i];
      return dialog;
    }
  }
}

export function setDialogInputValue(dialogTitle: string, inputName: string,
  value: string | number | boolean | DG.Column | DG.Column[]): void {
  returnDialog(dialogTitle)!.input(inputName).value = value;
}

export async function uploadProject(projectName: string, tableInfo: DG.TableInfo,
  view: DG.TableView, df: DG.DataFrame): Promise<void> {
  const project = DG.Project.create();
  project.name = projectName;
  project.addChild(tableInfo);
  project.addChild(view.saveLayout());  
  await grok.dapi.layouts.save(view.saveLayout());
  await grok.dapi.tables.uploadDataFrame(df);
  await grok.dapi.tables.save(tableInfo);
  await grok.dapi.projects.save(project);
}

export function findViewer(viewerName: string, view: DG.TableView): DG.Viewer | undefined {
  let viewer: DG.Viewer;
  for (let i: number = 0; i < Array.from(view.viewers).length; i++) {
    if (Array.from(view.viewers)[i].type == viewerName) {
      viewer = Array.from(view.viewers)[i];
      return viewer;
    }
  }
}

export function checkHTMLElementbyInnerText(className: string, innerText: string): void {
  const elements = document.getElementsByClassName(className);
  let check = false;
  let element;
  for (let i = 0; i < elements.length; i++ ) {
    element = elements[i] as HTMLElement;
    if (element.innerText == innerText)
      check = true;
  }
  if (check == false)
    throw new Error('element with innerText = "' + innerText + '" not found');
}

export function getHTMLElementbyInnerText(className: string, innerText: string): HTMLElement | undefined {
  const elements = document.getElementsByClassName(className);
  let element;
  for (let i = 0; i < elements.length; i++ ) {        
    element = elements[i] as HTMLElement;
    if (element.innerText == innerText)
      return element;
  }
}

export async function waitForElement(selector: string, error: string, wait=3000): Promise<HTMLElement> {
  return new Promise((resolve, reject) => {
    if (document.querySelector(selector))
      return resolve(document.querySelector(selector) as HTMLElement);

    const observer = new MutationObserver(() => {
      if (document.querySelector(selector)) {
        clearTimeout(timeout);
        observer.disconnect();
        resolve(document.querySelector(selector) as HTMLElement);
      }
    });
    
    const timeout = setTimeout(() => {
      observer.disconnect();
      reject(new Error(error));
    }, wait);

    observer.observe(document.body, {
      childList: true,
      subtree: true,
    });
  });
}
