import * as ui from 'datagrok-api/ui';

import '../../css/breadcrumbs.css';


export class Breadcrumbs {
  header: HTMLDivElement;
  content: HTMLDivElement;
  root: HTMLDivElement;
  activeTab: HTMLDivElement;
  disabledTabNames: string[];


  constructor(nameList: string[], contentList: HTMLElement[], disabledTabNames?: string[]) {
    this.header = document.createElement('div');
    this.content = document.createElement('div');
    this.root = document.createElement('div');
    this.activeTab = document.createElement('div');
    this.disabledTabNames = [];

    this._updateElements(nameList, contentList, disabledTabNames);
    this._initEventListeners();
  }


  private _updateElements(nameList: string[], contentList: HTMLElement[], disabledTabNames?: string[]) {
    const headerElements = nameList.map((name) => {
      const element = ui.div(name, `${name}-header-element breadcrumbs-header-element`);
      if (disabledTabNames?.includes(name))
        element.classList.add('disabled');
      return element;
    });
    this.header = ui.divH(headerElements, 'breadcrumbs-header-stripe');
    (this.header.firstChild as HTMLElement).classList.add('selected');
    this.activeTab = headerElements[0];
    this.disabledTabNames = disabledTabNames!;

    const contentElements = contentList.map((content, index) => ui.panel([content],
      `${nameList[index]}-content-element breadcrumbs-content-element`));
    this.content = ui.div(contentElements);
    (this.content.firstChild as HTMLElement).classList.add('visible');

    this.root = ui.div([this.header, this.content]);
  }

  private _initEventListeners() {
    const headerElements = this.header.getElementsByClassName('breadcrumbs-header-element');
    for (let i = 0; i < headerElements.length; i++) {
      const currentNode = (headerElements[i] as HTMLDivElement);
      const currentNodeName = (currentNode.firstChild as HTMLElement).innerHTML;

      this._addHeaderElementClickEvent(currentNode, currentNodeName);
    }
  }

  private _addHeaderElementClickEvent(element: HTMLDivElement, name: string) {
    element.addEventListener('click', () => {
      this.setActiveElement(name);
    }, false);
  }


  setActiveElement(name: string) {
    const currentHeaderElement = this.header.getElementsByClassName(`${name}-header-element`)![0];

    this.header.getElementsByClassName('selected')[0]?.classList.remove('selected');
    currentHeaderElement.classList.add('selected');

    this.content.getElementsByClassName('visible')[0]?.classList.remove('visible');
    this.content.getElementsByClassName(`${name}-content-element`)[0].classList.add('visible');

    this.activeTab = (currentHeaderElement as HTMLDivElement);
  }

  changeTabState(name: string, disabled: boolean) {
    const currentHeaderElement = this.header.getElementsByClassName(`${name}-header-element`)![0];
    if (disabled) {
      this.disabledTabNames[this.disabledTabNames.length] = name;
      currentHeaderElement.classList.add('disabled');
      return;
    }

    const index = this.disabledTabNames.indexOf(name);
    this.disabledTabNames.splice(index, 1);
    currentHeaderElement.classList.remove('disabled');
  }

  getTabContent(name: string): HTMLElement {
    return (this.content.getElementsByClassName(`${name}-content-element`)[0] as HTMLElement);
  }

  editTab(tabName: string, newName: string, content: HTMLElement) {
    if (this.disabledTabNames.includes(tabName)) {
      const index = this.disabledTabNames.indexOf(tabName);
      this.disabledTabNames.splice(index, 1);
    }

    const currentHeaderElement = this.header.getElementsByClassName(`${tabName}-header-element`)![0];
    const newHeaderElement = ui.div(newName, `${newName}-header-element breadcrumbs-header-element`);
    this.header.replaceChild(newHeaderElement, currentHeaderElement);

    const currentContentElement = this.content.getElementsByClassName(`${tabName}-content-element`)![0];
    const newContentElement = ui.panel([content], `${newName}-content-element breadcrumbs-content-element`);
    this.content.replaceChild(newContentElement, currentContentElement);

    this._addHeaderElementClickEvent(newHeaderElement, newName);
  }

  appendTab(name: string, content: HTMLElement, disabled?: boolean) {
    const newHeaderElement = ui.div(name, `${name}-header-element breadcrumbs-header-element`);
    this.header.appendChild(newHeaderElement);
    if (disabled) {
      this.disabledTabNames[this.disabledTabNames.length] = name;
      newHeaderElement.classList.add('disabled');
    }

    const newContentElement = ui.panel([content], `${name}-content-element breadcrumbs-content-element`);
    this.content.appendChild(newContentElement);

    this._addHeaderElementClickEvent(newHeaderElement, name);
  }


  onTabChanged(callback: Function) {
    const headerElements = this.header.getElementsByClassName('breadcrumbs-header-element');
    for (let i = 0; i < headerElements.length; i++) {
      const currentNode = headerElements[i];
      currentNode.addEventListener('click', () => {
        callback();
      }, false);
    }
  }
}
