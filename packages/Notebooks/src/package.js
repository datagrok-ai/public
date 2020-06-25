import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from "datagrok-api/dg";

export let _package = new DG.Package();

import '@jupyterlab/application/style/index.css';
import '@jupyterlab/codemirror/style/index.css';
import '@jupyterlab/completer/style/index.css';
import '@jupyterlab/documentsearch/style/index.css';
import '@jupyterlab/notebook/style/index.css';
import '../css/application-base.css'
import '../css/notebooks.css';
import '../css/theme-light-extension-index.css';
import '../css/ui-components-base.css';

import { PageConfig } from '@jupyterlab/coreutils';
import { CommandRegistry } from '@lumino/commands';
import { Widget} from '@lumino/widgets';
import { ServiceManager, ServerConnection } from '@jupyterlab/services';
import { MathJaxTypesetter } from '@jupyterlab/mathjax2';
import { NotebookPanel, NotebookWidgetFactory, NotebookModelFactory, NotebookActions, CellTypeSwitcher } from '@jupyterlab/notebook';
import { CompleterModel, Completer, CompletionHandler, KernelConnector } from '@jupyterlab/completer';
import { editorServices } from '@jupyterlab/codemirror';
import { DocumentManager } from '@jupyterlab/docmanager';
import { DocumentRegistry } from '@jupyterlab/docregistry';
import { RenderMimeRegistry, standardRendererFactories as initialFactories } from '@jupyterlab/rendermime';
import { SetupCommands } from './commands';
import {sessionContextDialogs} from "@jupyterlab/apputils";


class NotebookView extends DG.ViewBase {
    constructor(params, path) {
        super(params, path);

        this.TYPE = 'Notebook';
        this.PATH = '/notebook';

        this.openAsHtmlIcon = ui.bigButton('HTML', () => { this.htmlMode().then() }, 'Open as HTML', true);
        this.openAsHtmlIcon.style.backgroundColor = 'var(--green-2)';
        this.editIcon = ui.bigButton('EDIT', () => { this.editMode().then() }, 'Edit notebook', true);
        this.editIcon.style.backgroundColor = 'var(--green-2)';

        this.html = '';
        this.saveAsComboPopup = ui.comboPopup(ui.iconFA('arrow-to-bottom', () => {}), ['As HTML', 'As PDF'], (item) => {
            if (item === 'As HTML') {
                let a = document.createElement('a');
                a.download = `${this.notebook.name}.html`;
                a.href = URL.createObjectURL(new Blob([this.html], {type: 'text/plain'})).toString();
                a.click();
            } else {
                let w = window.open('', '','height=650,width=900,top=100,left=150');
                let d = w['document'];
                d['body']['innerHTML'] = this.html;
                d['title'] = this.notebook.name;
                w.print();
            }
        });

        this.environmentInput = ui.choiceInput('Environment:', 'default', ['default', 'python3']);

        this.id = ('id' in params) ? params['id']
            : ((path !== null && path !== undefined) ? path.replace(`${this.PATH}/`, '')
            : null);

        if ('file' in params)
            this.notebookFile = params['file'];

        let edit = 'edit' in params ? params['edit'] : false;
        let html = 'html' in params ? params['html'] : null;

        if (edit)
            this.editMode().then();
        else
            this.htmlMode(html).then();
    }

    get type() { return this.TYPE };
    get helpUrl() { return '/help/compute/jupyter-notebook.md'; }
    get name() { return (this.notebook !== null && this.notebook !== undefined) ? this.notebook.name : 'Notebook' };
    get path() { return `${this.PATH}/${this.id}` };

    getIcon() {
        let img = document.createElement('img');
        img.src = '/images/entities/jupyter.png';
        img.height = 18;
        img.width = 18;
        return img;
    };

    saveStateMap() { return {'notebookId': this.id }; }
    loadStateMap(stateMap) { open(stateMap['notebookId']); }

    handlePath(path) {
        let id = path.replace(`${this.PATH}/`, '');
        open(id);
    }

    acceptsPath(path) { return path.startsWith(this.PATH); }

    async initNotebook() {
        this.notebook = await grok.dapi.notebooks.find(this.id);
    }

    async htmlMode(html = null) {
        this.html = html;
        await this.initNotebook();
        removeChildren(this.root);
        if (this.html === null)
            this.html = await this.notebook.toHtml();
        let iframe = document.createElement('iframe');
        iframe.src = 'data:text/html;base64,' + btoa(this.html);
        iframe.classList.add("grok-jupyter-notebook-html-iframe");
        let container = ui.div([iframe], 'grok-jupyter-notebook-container');
        let view = ui.div([container], 'd4-root,d4-flex-col');
        this.setRibbonPanels([[this.saveAsComboPopup, this.editIcon]], true);
        this.root.appendChild(view);
    }

    async editMode() {
        await this.initNotebook();
        removeChildren(this.root);

        let notebookPath = this.notebookFile;
        const manager = new ServiceManager({serverSettings: NotebookView.getSettings()});
        await manager.ready;

        // Initialize the command registry with the bindings.
        // Setup the keydown listener for the document.
        const commands = new CommandRegistry();
        const useCapture = true;
        document.addEventListener('keydown', event => { commands.processKeydownEvent(event); }, useCapture);

        const renderMime = new RenderMimeRegistry({
            initialFactories: initialFactories,
            latexTypesetter: new MathJaxTypesetter({
                url: PageConfig.getOption('mathjaxUrl'),
                config: PageConfig.getOption('mathjaxConfig')
            })
        });

        const opener = { open: (widget) => {}};
        const docRegistry = new DocumentRegistry();
        const docManager = new DocumentManager({registry: docRegistry, manager, opener});
        const mFactory = new NotebookModelFactory({});
        const editorFactory = editorServices.factoryService.newInlineEditor;
        const contentFactory = new NotebookPanel.ContentFactory({editorFactory});

        const wFactory = new NotebookWidgetFactory({
            name: 'Notebook',
            modelName: 'notebook',
            fileTypes: ['notebook'],
            defaultFor: ['notebook'],
            preferKernel: true,
            canStartKernel: true,
            rendermime: renderMime,
            contentFactory,
            mimeTypeService: editorServices.mimeTypeService
        });
        docRegistry.addModelFactory(mFactory);
        docRegistry.addWidgetFactory(wFactory);

        const nbWidget = docManager.open(notebookPath);

        const editor = nbWidget.content.activeCell && nbWidget.content.activeCell.editor;
        const model = new CompleterModel();
        const completer = new Completer({editor, model});
        const sessionContext = nbWidget.context.sessionContext;
        const connector = new KernelConnector({session: sessionContext.session});
        const handler = new CompletionHandler({completer, connector});

        void sessionContext.ready.then(() => {
            handler.connector = new KernelConnector({session: sessionContext.session});
        });

        handler.editor = editor;

        nbWidget.content.activeCellChanged.connect((sender, cell) => {
            handler.editor = cell !== null && cell !== undefined ? cell.editor : null;
        });

        completer.hide();

        Widget.attach(nbWidget, this.root);
        Widget.attach(completer, this.root);

        this.setRibbonPanels([
            [ this.openAsHtmlIcon ],
            [
                ui.iconFA('save', () => nbWidget.context.save(), 'Save notebook'),
                ui.iconFA('plus', () => NotebookActions.insertBelow(nbWidget.content), 'Insert a cell before'),
                ui.iconFA('cut', () => NotebookActions.cut(nbWidget.content), 'Cut cell'),
                ui.iconFA('copy', () => NotebookActions.copy(nbWidget.content), 'Copy cell'),
                ui.iconFA('paste', () => NotebookActions.copy(nbWidget.content), 'Paste cell'),
                ui.iconFA('play', () => NotebookActions.runAndAdvance(nbWidget.content, nbWidget.context.sessionContext), 'Run cell'),
                ui.iconFA('stop', () => nbWidget.context.sessionContext.session.kernel.interrupt(), 'Interrupt Kernel'),
                ui.iconFA('redo', () => sessionContextDialogs.restart(nbWidget.context.sessionContext), 'Restart Kernel'),
                ui.iconFA('forward', () => sessionContextDialogs.restart(nbWidget.context.sessionContext)
                    .then(restarted => { if (restarted) NotebookActions.runAll(nbWidget.content, nbWidget.context.sessionContext); }),
                    'Restart Kernel and run all cells'),
                new CellTypeSwitcher(nbWidget.content).node,
            ],
            [ this.environmentInput.root ]
        ], true);
        nbWidget.toolbar.hide();

        //ui.tools.handleResize(this.root, (w, h) => nbWidget.update());
        this.root.classList.add('grok-notebook-view');

        SetupCommands(commands, nbWidget, handler);
    }

    static getSettings() {
        const _settings = {
            baseUrl: grok.settings.jupyterNotebook,
            token: grok.settings.jupyterNotebookToken,
            mathjaxUrl: 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js',
            mathjaxConfig: 'TeX-AMS_CHTML-full,Safe'
        };

        for (let key in _settings) PageConfig.setOption(key, _settings[key]);

        let settings = ServerConnection.defaultSettings;
        settings.baseUrl = PageConfig.getBaseUrl();
        settings.wsUrl = PageConfig.getWsUrl();
        settings.token = PageConfig.getToken();

        return settings;
    }
}


//name: Notebook
//description: Creates a Notebook View
//input: map params =
//input: string path =
//tags: view
//output: view result
export function notebookView(params = null, path = '') {
    return new NotebookView(params, path);
}


function removeChildren(node) {
    while (node.firstChild)
        node.removeChild(node.firstChild);
}
