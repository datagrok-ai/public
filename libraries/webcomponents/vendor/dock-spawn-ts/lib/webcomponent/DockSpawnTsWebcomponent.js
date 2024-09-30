import { DockManager } from '../DockManager.js';
import { PanelContainer } from '../PanelContainer.js';
import { PanelType } from '../enums/PanelType.js';
import { style, style1, style2 } from './styles';
export class DockSpawnTsWebcomponent extends HTMLElement {
    dockManager;
    slotId = 0;
    windowResizedBound;
    slotElementMap;
    observer;
    resizeObserver;
    initialized = false;
    initFinished = false;
    elementContainerMap = new Map();
    constructor() {
        super();
        const shadowRoot = this.attachShadow({ mode: 'open' });
        shadowRoot.adoptedStyleSheets = [style, style1, style2];
        this.windowResizedBound = this.windowResized.bind(this);
        this.slotElementMap = new WeakMap();
    }
    saveLayout() {
        return this.dockManager.saveState();
    }
    async loadLayout(serializedState) {
        this.observer.disconnect();
        const oldPanels = this.dockManager.getPanels();
        await this.dockManager.loadState(serializedState);
        this.observer.observe(this, { childList: true });
        oldPanels.forEach((panel) => panel.elementContentContainer.remove());
        this.dockManager.element.firstChild.remove();
        this.dockManager.invalidate();
    }
    initDockspawn() {
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
                    const slot = dockNode.elementContent;
                    const element = this.slotElementMap.get(slot);
                    if (element) {
                        if (this.contains(element))
                            this.removeChild(element);
                        this.dispatchEvent(new CustomEvent('panel-closed', { detail: element }));
                    }
                },
            });
            for (const element of this.children)
                this.handleAddedChildNode(element);
            this.observer = new MutationObserver((mutations) => {
                mutations.forEach((mutation) => {
                    mutation.addedNodes.forEach((node) => {
                        this.handleAddedChildNode(node);
                    });
                    mutation.removedNodes.forEach((node) => {
                        this.handleRemovedChildNode(node);
                    });
                });
            });
            this.observer.observe(this, { childList: true });
            this.resizeObserver = new ResizeObserver(() => {
                this.resize();
            });
            this.resizeObserver.observe(this);
            this.initFinished = true;
        }, 50);
    }
    getElementInSlot(slot) {
        return this.slotElementMap.get(slot);
    }
    handleAddedChildNode(element) {
        if (element instanceof Comment)
            return;
        const slot = document.createElement('slot');
        const slotName = 'slot_' + this.slotId++;
        slot.name = slotName;
        let dockPanelType = PanelType.panel;
        const dockPanelTypeAttribute = element.getAttribute('dock-spawn-panel-type');
        if (dockPanelTypeAttribute)
            dockPanelType = dockPanelTypeAttribute;
        const hideCloseButton = element.hasAttribute('dock-spawn-hide-close-button');
        const title = element.getAttribute('dock-spawn-title');
        const container = new PanelContainer(slot, this.dockManager, title, dockPanelType, hideCloseButton);
        element.slot = slotName;
        this.slotElementMap.set(slot, element);
        this.elementContainerMap.set(element, container);
        let dockRatio = 0.5;
        const dockRatioAttribute = element.getAttribute('dock-spawn-dock-ratio');
        if (dockRatioAttribute)
            dockRatio = dockRatioAttribute;
        const dockType = element.getAttribute('dock-spawn-dock-type');
        let dockRelativeTo = this.dockManager.context.model.documentManagerNode;
        const dockToAttribute = element.getAttribute('dock-spawn-dock-to');
        if (dockToAttribute) {
            //@ts-ignore
            const dockToElement = this.getRootNode().getElementById(dockToAttribute);
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
        if (element.style.display == 'none')
            element.style.display = 'block';
    }
    handleRemovedChildNode(element) {
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
        if (this.initFinished)
            this.resize();
    }
    disconnectedCallback() {
        window.removeEventListener('resize', this.windowResizedBound);
        window.removeEventListener('orientationchange', this.windowResizedBound);
    }
    windowResized() {
        this.resize();
    }
    resize() {
        this.dockManager.resize(this.clientWidth, this.clientHeight);
    }
    getDockNodeForElement(elementOrContainer) {
        let element = elementOrContainer;
        if (element.containerElement)
            element = elementOrContainer.containerElement;
        return this.dockManager.findNodeFromContainerElement(element);
    }
    dockFill(element, panelType, dockNode, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockFill(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container);
    }
    dockLeft(element, panelType, dockNode, ratio, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockLeft(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }
    dockRight(element, panelType, dockNode, ratio, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockRight(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }
    dockUp(element, panelType, dockNode, ratio, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockUp(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }
    dockDown(element, panelType, dockNode, ratio, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockDown(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }
    floatDialog(element, x, y, width, height, panelType, title, hideCloseButton) {
        const container = new PanelContainer(element, this.dockManager, title, panelType, hideCloseButton);
        const dlg = this.dockManager.floatDialog(container, x, y, null);
        dlg.resize(width, height);
    }
}
