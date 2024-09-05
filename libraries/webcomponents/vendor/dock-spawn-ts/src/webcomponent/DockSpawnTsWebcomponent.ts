import {DockManager} from '../DockManager.js';
import {PanelContainer} from '../PanelContainer.js';
import {PanelType} from '../enums/PanelType.js';
import {DockNode} from '../DockNode.js';
import {style, style1, style2} from './styles';

export class DockSpawnTsWebcomponent extends HTMLElement {
  public dockManager: DockManager;
  private slotId: number = 0;
  private windowResizedBound;
  private slotElementMap: WeakMap<HTMLSlotElement, HTMLElement>;
  private observer: MutationObserver;
  private initialized = false;
  private elementContainerMap: Map<HTMLElement, PanelContainer> = new Map();

  constructor() {
    super();

    const shadowRoot = this.attachShadow({mode: 'open'});
    shadowRoot.adoptedStyleSheets = [style, style1, style2];

    this.windowResizedBound = this.windowResized.bind(this);
    this.slotElementMap = new WeakMap();
  }

  private initDockspawn() {
    const dockSpawnDiv = document.createElement('div');
    dockSpawnDiv.id = 'dockSpawnDiv';
    dockSpawnDiv.style.width = '100%';
    dockSpawnDiv.style.height = '100%';
    dockSpawnDiv.style.position = 'relative';
    this.shadowRoot.appendChild(dockSpawnDiv);

    this.dockManager = new DockManager(dockSpawnDiv);
    this.dockManager.config.dialogRootElement = dockSpawnDiv;

    setTimeout(() => {
      this.dockManager.initialize();

      this.dockManager.addLayoutListener({
        onClosePanel: (dockManager, dockNode) => {
          const slot = dockNode.elementContent as any as HTMLSlotElement;
          const element = this.slotElementMap.get(slot);
          if (element) {
            if (this.contains(element)) this.removeChild(element);
            this.dispatchEvent(new CustomEvent('panel-closed', {detail: element}));
          }
        },
      });

      this.onresize = () => {
        this.dockManager.resize(this.clientWidth, this.clientHeight);
      };

      for (const element of this.children)
        this.handleAddedChildNode(element as HTMLElement);


      this.observer = new MutationObserver((mutations) => {
        mutations.forEach((mutation) => {
          mutation.addedNodes.forEach((node) => {
            this.handleAddedChildNode(node as HTMLElement);
          });
          mutation.removedNodes.forEach((node) => {
            this.handleRemovedChildNode(node as HTMLElement);
          });
        });
      });
      this.observer.observe(this, {childList: true});
    }, 50);
  }

  public getElementInSlot(slot: HTMLSlotElement): HTMLElement {
    return this.slotElementMap.get(slot);
  }

  private handleAddedChildNode(element: HTMLElement) {
    if (element instanceof Comment) return;

    const slot = document.createElement('slot');
    const slotName = 'slot_' + this.slotId++;
    slot.name = slotName;

    let dockPanelType = PanelType.panel;
    const dockPanelTypeAttribute = element.getAttribute('dock-spawn-panel-type');
    if (dockPanelTypeAttribute)
      dockPanelType = <PanelType><any>dockPanelTypeAttribute;
    const hideCloseButton = element.hasAttribute('dock-spawn-hide-close-button');
    const container = new PanelContainer(slot, this.dockManager, element.title, dockPanelType, hideCloseButton);
    element.slot = slotName;
    this.slotElementMap.set(slot, (<HTMLElement>element));
    this.elementContainerMap.set(element, container);

    let dockRatio: number = 0.5;
    const dockRatioAttribute = element.getAttribute('dock-spawn-dock-ratio');
    if (dockRatioAttribute)
      dockRatio = <number><any>dockRatioAttribute;
    const dockType = element.getAttribute('dock-spawn-dock-type');

    let dockRelativeTo = this.dockManager.context.model.documentManagerNode;
    const dockToAttribute = element.getAttribute('dock-spawn-dock-to');
    if (dockToAttribute) {
      //@ts-ignore
      const dockToElement = this.getRootNode().getElementById(dockToAttribute) as HTMLElement;
      dockRelativeTo = this.dockManager.findNodeFromContainerElement(this.elementContainerMap.get(dockToElement).containerElement);
    }

    if (dockType == 'left')
      this.dockManager.dockLeft(dockRelativeTo, container, dockRatio);
    else if (dockType == 'right')
      this.dockManager.dockRight(dockRelativeTo, container, dockRatio);
    else if (dockType == 'up')
      this.dockManager.dockUp(dockRelativeTo, container, dockRatio);
    else if (dockType == 'down')
      this.dockManager.dockDown(dockRelativeTo, container, dockRatio);
    else
      this.dockManager.dockFill(dockRelativeTo, container);

    if ((<HTMLElement>element).style.display == 'none')
      (<HTMLElement>element).style.display = 'block';
  }

  private handleRemovedChildNode(element: HTMLElement) {
    const node = this.elementContainerMap.get(element);
    if (node)
      node.close();
    this.elementContainerMap.delete(element);
  }

  connectedCallback() {
    if (!this.initialized) {
      this.initDockspawn();
      this.initialized = true;
    }
    window.addEventListener('resize', this.windowResizedBound);
    window.addEventListener('orientationchange', this.windowResizedBound);
  }

  disconnectedCallback() {
    window.removeEventListener('resize', this.windowResizedBound);
    window.removeEventListener('orientationchange', this.windowResizedBound);
  }

  private windowResized() {
    this.resize();
  }

  resize() {
    this.dockManager.resize(this.clientWidth, this.clientHeight);
  }

  getDockNodeForElement(elementOrContainer: HTMLElement | PanelContainer): DockNode {
    let element = elementOrContainer as HTMLElement;
    if ((<any>element).containerElement)
      element = (<any>elementOrContainer).containerElement as HTMLElement;
    return this.dockManager.findNodeFromContainerElement(element);
  }

  dockFill(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    this.dockManager.dockFill(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container);
  }

  dockLeft(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    this.dockManager.dockLeft(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
  }

  dockRight(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    this.dockManager.dockRight(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
  }

  dockUp(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    this.dockManager.dockUp(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
  }

  dockDown(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    this.dockManager.dockDown(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
  }

  floatDialog(element: HTMLElement, x: number, y: number, width: number, height: number, panelType?: PanelType, title?: string, hideCloseButton?: boolean) {
    const container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
    const dlg = this.dockManager.floatDialog(container, x, y, null);
    dlg.resize(width, height);
  }
}

window.customElements.define('dock-spawn-ts', DockSpawnTsWebcomponent);
