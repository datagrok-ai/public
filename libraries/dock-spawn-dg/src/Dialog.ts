import { DockManager } from "./DockManager";
import { Point } from "./Point";
import { PanelContainer } from "./PanelContainer";
import { DraggableContainer } from "./DraggableContainer";
import { ResizableContainer } from "./ResizableContainer";
import { EventHandler } from "./EventHandler";
import { Utils } from "./Utils";
import { Localizer } from "./i18n/Localizer";
import { DockNode } from "./DockNode";

export class Dialog {
    elementDialog: HTMLDivElement & { floatingDialog: Dialog };
    draggable: DraggableContainer;
    panel: PanelContainer;
    dockManager: DockManager;
    eventListener: DockManager;
    position: Point;
    resizable: ResizableContainer;
    disableResize: boolean;
    mouseDownHandler: any;
    onKeyPressBound: any;
    noDocking: boolean;
    isHidden: boolean;
    keyPressHandler: EventHandler;
    focusHandler: EventHandler;
    grayoutParent: PanelContainer;
    contextmenuHandler: EventHandler;
    _ctxMenu: HTMLDivElement;
    _windowsContextMenuCloseBound: any;

    constructor(panel: PanelContainer, dockManager: DockManager, grayoutParent?: PanelContainer, disableResize?: boolean) {
        this.panel = panel;
        this.dockManager = dockManager;
        this.eventListener = dockManager;
        this.grayoutParent = grayoutParent;
        this.disableResize = disableResize;
        this._initialize();
        this.dockManager.context.model.dialogs.push(this);
        this.position = dockManager.defaultDialogPosition;
        this.dockManager.notifyOnCreateDialog(this);
        panel.isDialog = true;
    }

    saveState(x: number, y: number) {
        this.position = new Point(x, y);
        this.dockManager.notifyOnChangeDialogPosition(this, x, y);
    }

    static fromElement(id: string, dockManager: DockManager) {
        return new Dialog(new PanelContainer(<HTMLElement>document.getElementById(id), dockManager), dockManager, null);
    }

    _initialize() {
        this.panel.floatingDialog = this;
        this.elementDialog = Object.assign(document.createElement('div'), { floatingDialog: this });
        this.elementDialog.tabIndex = 0;
        this.elementDialog.appendChild(this.panel.elementPanel);
        this.draggable = new DraggableContainer(this, this.panel, this.elementDialog, this.panel.elementTitle);
        this.resizable = new ResizableContainer(this, this.draggable, this.draggable.topLevelElement, this.disableResize);

        this.dockManager.config.dialogRootElement.appendChild(this.elementDialog);
        this.elementDialog.classList.add('dialog-floating');

        this.focusHandler = new EventHandler(this.elementDialog, 'focus', this.onFocus.bind(this), true);
        this.mouseDownHandler = new EventHandler(this.elementDialog, 'pointerdown', this.onMouseDown.bind(this), true);
        this.keyPressHandler = new EventHandler(this.elementDialog, 'keypress', this.dockManager.onKeyPressBound, true);
        this.contextmenuHandler = new EventHandler(this.panel.elementTitle, 'contextmenu', this.oncontextMenuClicked.bind(this));

        this.resize(this.panel.elementPanel.clientWidth, this.panel.elementPanel.clientHeight);
        this.isHidden = false;

        if (this.grayoutParent != null) {
            this.grayoutParent.grayOut(true);
        }
        this.bringToFront();
    }

    setPosition(x: number, y: number) {
        let rect = this.dockManager.config.dialogRootElement.getBoundingClientRect();
        this.position = new Point(x - rect.left, y - rect.top);
        this.elementDialog.style.left = (x - rect.left) + 'px';
        this.elementDialog.style.top = (y - rect.top) + 'px';
        this.panel.setDialogPosition(x, y);
        this.dockManager.notifyOnChangeDialogPosition(this, x, y);
    }

    getPosition(): Point {
        return new Point(this.position ? this.position.x : 0, this.position ? this.position.y : 0);
    }

    onFocus() {
        if (this.dockManager.activePanel != this.panel)
            this.dockManager.activePanel = this.panel;
    }

    onMouseDown(e: PointerEvent) {
        if (e.button != 2)
            this.bringToFront();
    }

    destroy() {
        this.panel.lastDialogSize = { width: this.resizable.width, height: this.resizable.height };

        if (this.focusHandler) {
            this.focusHandler.cancel();
            delete this.focusHandler;
        }
        if (this.mouseDownHandler) {
            this.mouseDownHandler.cancel();
            delete this.mouseDownHandler;
        }
        if (this.keyPressHandler) {
            this.keyPressHandler.cancel();
            delete this.keyPressHandler;
        }
        if (this.contextmenuHandler) {
            this.contextmenuHandler.cancel();
            delete this.contextmenuHandler;
        }
        Utils.removeNode(this.elementDialog);
        this.draggable.removeDecorator();
        Utils.removeNode(this.panel.elementPanel);
        Utils.arrayRemove(this.dockManager.context.model.dialogs, this);
        delete this.panel.floatingDialog;

        if (this.grayoutParent) {
            this.grayoutParent.grayOut(false);
        }
    }

    resize(width: number, height: number) {
        this.resizable.resize(width, height);
    }

    setTitle(title: string) {
        this.panel.setTitle(title);
    }

    setTitleIcon(iconName: string) {
        this.panel.setTitleIcon(iconName);
    }

    bringToFront() {
        this.panel.elementContentContainer.style.zIndex = <any>this.dockManager.zIndexDialogCounter++;
        this.elementDialog.style.zIndex = <any>this.dockManager.zIndexDialogCounter++;
        this.dockManager.activePanel = this.panel;
    }

    hide() {
        this.elementDialog.style.zIndex = '0';
        this.panel.elementContentContainer.style.zIndex = '';
        this.elementDialog.style.display = 'none';
        if (!this.isHidden) {
            this.isHidden = true;
            this.dockManager.notifyOnHideDialog(this);
        }
        if (this.grayoutParent) {
            this.grayoutParent.grayOut(false);
        }
    }

    close() {
        this.hide();
        this.remove();
        this.dockManager.notifyOnClosePanel(this.panel);
        this.destroy();
    }

    remove() {
        this.elementDialog.parentNode.removeChild(this.elementDialog);
    }

    show() {
        this.panel.elementContentContainer.style.zIndex = <any>this.dockManager.zIndexDialogCounter++;
        this.elementDialog.style.zIndex = <any>this.dockManager.zIndexDialogCounter++;
        this.elementDialog.style.display = 'block';
        if (this.isHidden) {
            this.isHidden = false;
            this.dockManager.notifyOnShowDialog(this);
        }
    }

    static createContextMenuContentCallback = (dialog: Dialog, documentMangerNodes: DockNode[]): Node[] => {
        if (!dialog.panel._hideCloseButton) {
            return [];
        }

        const result = [];

        let btnCloseDialog = document.createElement('div');
        btnCloseDialog.innerText = Localizer.getString('CloseDialog');
        result.push(btnCloseDialog);

        btnCloseDialog.onclick = () => {
            dialog.panel.close();
            dialog.closeContextMenu();
        };

        if (dialog.dockManager.config.enableBrowserWindows) {
            let btnNewBrowserWindow = document.createElement('div');
            btnNewBrowserWindow.innerText = Localizer.getString('NewBrowserWindow');
            result.push(btnNewBrowserWindow);

            btnNewBrowserWindow.onclick = () => {
                dialog.panel.undockToBrowserDialog();
                dialog.closeContextMenu();
            };
        }

        return result;
    }

    oncontextMenuClicked(e: MouseEvent) {
        e.preventDefault();

        if (!this._ctxMenu && Dialog.createContextMenuContentCallback) {
            const menuItems = Dialog.createContextMenuContentCallback(
                this,
                this.dockManager.context.model.documentManagerNode.children
            );

            if (menuItems.length == 0) {
                return;
            }

            this._ctxMenu = document.createElement('div');
            this._ctxMenu.className = 'dockspab-tab-handle-context-menu';
            this._ctxMenu.append(...menuItems);
            this._ctxMenu.style.left = e.pageX + "px";
            this._ctxMenu.style.top = e.pageY + "px";
            document.body.appendChild(this._ctxMenu);
            this._windowsContextMenuCloseBound = this.windowsContextMenuClose.bind(this)
            window.addEventListener('pointerup', this._windowsContextMenuCloseBound);
        } else {
            this.closeContextMenu();
        }
    }

    closeContextMenu() {
        if (this._ctxMenu) {
            document.body.removeChild(this._ctxMenu);
            delete this._ctxMenu;
            window.removeEventListener('pointerup', this._windowsContextMenuCloseBound);
        }
    }

    windowsContextMenuClose(e: Event) {
        let cp = e.composedPath();
        for (let i in cp) {
            let el = cp[i];
            if (el == this._ctxMenu)
                return;
        }
        this.closeContextMenu();
    }
}
