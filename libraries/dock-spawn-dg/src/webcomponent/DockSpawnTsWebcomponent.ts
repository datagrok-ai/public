import { DockManager } from "../DockManager";
import { PanelContainer } from "../PanelContainer";
import { PanelType } from "../enums/PanelType";
import { DockNode } from "../DockNode";
import { faStyle, style, style1, style2 } from "./styles";
import { Observable, Subject, Subscriber } from "rxjs";
import { debounceTime, takeUntil } from "rxjs/operators";

const elementResized = (elem: Element) => {
    return new Observable((subscriber: Subscriber<ResizeObserverEntry[]>) => {
        const resizeObserver = new ResizeObserver((entries: ResizeObserverEntry[]) => subscriber.next(entries),);
        resizeObserver.observe(elem);
        return () => resizeObserver.disconnect();
    });
};

export class DockSpawnTsWebcomponent extends HTMLElement {
    public dockManager: DockManager;
    private slotId: number = 0;
    private slotElementMap: WeakMap<HTMLSlotElement, HTMLElement>;
    private observer: MutationObserver;
    private initialized = false;
    private elementContainerMap: Map<HTMLElement, PanelContainer> = new Map();
    private destroy$ = new Subject<true>();

    constructor() {
        super();

        const shadowRoot = this.attachShadow({ mode: 'open' });
        shadowRoot.adoptedStyleSheets = [style, style1, style2, faStyle];

        this.slotElementMap = new WeakMap();
    }

    private initDockspawn() {
        const dockSpawnDiv = document.createElement('div')
        dockSpawnDiv.id = "dockSpawnDiv";
        dockSpawnDiv.style.width = "100%";
        dockSpawnDiv.style.height = "100%";
        dockSpawnDiv.style.position = "relative";
        this.shadowRoot.appendChild(dockSpawnDiv);

        this.dockManager = new DockManager(dockSpawnDiv);
        this.dockManager.config.dialogRootElement = dockSpawnDiv;
    }

    private initDockManager() {
        this.dockManager.initialize();

        this.dockManager.addLayoutListener({
            onActivePanelChange: (_, panel, prevPanel) => {
                this.dispatchEvent(new CustomEvent('active-panel-changed', {detail: {
                    newPanel: panel?.title,
                    prevPanel: prevPanel?.title
                }}));
            },
            onClosePanel: (dockManager, dockNode) => {
                let slot = dockNode.elementContent as any as HTMLSlotElement;
                let element = this.slotElementMap.get(slot);
                if (element) {
                    if (this.contains(element)) {
                        this.observer.disconnect();
                        this.removeChild(element);
                        this.observer.observe(this, {childList: true});
                    }
                    this.dispatchEvent(new CustomEvent('panel-closed', {detail: element}));
                }
            },
        });

        for (const element of this.children) {
            this.handleAddedChildNode(element as HTMLElement);
        }

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

        elementResized(this).pipe(
            debounceTime(50),
            takeUntil(this.destroy$)
        ).subscribe(() => this.resize());
        this.dispatchEvent(new CustomEvent('manager-init-finished'));
    }

    public set activePanelTitle(panelTitle: string) {
        const foundPanel = this.dockManager?.getPanels().find((panel) => panel.title === panelTitle);
        if (foundPanel)
            this.dockManager.activePanel = foundPanel;
    }

    public getElementInSlot(slot: HTMLSlotElement): HTMLElement {
        return this.slotElementMap.get(slot);
    }

    private handleAddedChildNode(element: HTMLElement) {
        if (element instanceof Comment || (element instanceof Text && element.textContent.length === 0)) return;

        let slot = document.createElement('slot');
        let slotName = 'slot_' + this.slotId++;
        slot.name = slotName;

        let dockPanelType = PanelType.panel;
        let dockPanelTypeAttribute = element.getAttribute('dock-spawn-panel-type');
        if (dockPanelTypeAttribute)
            dockPanelType = <PanelType><any>dockPanelTypeAttribute;
        let hideCloseButton = element.hasAttribute('dock-spawn-hide-close-button');
        let title = element.getAttribute('dock-spawn-title');
        let panelIcon = element.getAttribute('dock-spawn-panel-icon');
        let zIndex = element.getAttribute('dock-spawn-z-index');
        let container = new PanelContainer(slot, this.dockManager, title, dockPanelType, hideCloseButton, panelIcon, zIndex);
        element.slot = slotName;
        this.slotElementMap.set(slot, (<HTMLElement>element));
        this.elementContainerMap.set(element, container);

        let dockRatio: number = 0.5;
        let dockRatioAttribute = element.getAttribute('dock-spawn-dock-ratio');
        if (dockRatioAttribute)
            dockRatio = <number><any>dockRatioAttribute;
        let dockType = element.getAttribute('dock-spawn-dock-type');

        let dockRelativeTo = this.dockManager.context.model.documentManagerNode;
        let dockToAttribute = element.getAttribute('dock-spawn-dock-to');
        if (dockToAttribute) {
            let dockToElement = (this.getRootNode() as HTMLElement)
                .querySelector(`[dock-spawn-title="${dockToAttribute}"]`) as HTMLElement | null;
            if (dockToElement && this.elementContainerMap.get(dockToElement)) {
                dockRelativeTo = this.dockManager
                    .findNodeFromContainerElement(this.elementContainerMap.get(dockToElement).containerElement) ?? dockRelativeTo;
            }
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
        let panel = this.elementContainerMap.get(element);
        if (panel)
            panel.close();
        this.elementContainerMap.delete(element);
    }

    connectedCallback() {
        if (!this.initialized) {
            this.initialized = true;
            this.initDockspawn();
            setTimeout(() => this.initDockManager(), 0);
            window.addEventListener('orientationchange', () => this.resize());
        }
    }

    disconnectedCallback() {}

    resize() {
        if (this.clientWidth > 0 && this.clientHeight > 0)
            this.dockManager.resize(this.clientWidth, this.clientHeight);
    }

    getDockNodeForElement(elementOrContainer: HTMLElement | PanelContainer): DockNode {
        let element = elementOrContainer as HTMLElement;
        if ((<any>element).containerElement)
            element = (<any>elementOrContainer).containerElement as HTMLElement;
        return this.dockManager.findNodeFromContainerElement(element);
    }

    dockFill(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockFill(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container);
    }

    dockLeft(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockLeft(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }

    dockRight(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockRight(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }

    dockUp(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockUp(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }

    dockDown(element: HTMLElement, panelType?: PanelType, dockNode?: DockNode, ratio?: number, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        this.dockManager.dockDown(dockNode != null ? dockNode : this.dockManager.context.model.documentManagerNode, container, ratio);
    }

    floatDialog(element: HTMLElement, x: number, y: number, width: number, height: number, panelType?: PanelType, title?: string, hideCloseButton?: boolean) {
        let container = new PanelContainer(element as HTMLElement, this.dockManager, title, panelType, hideCloseButton);
        let dlg = this.dockManager.floatDialog(container, x, y, null);
        dlg.resize(width, height);
    }
}
