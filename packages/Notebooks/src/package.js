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
import { CommandPalette, SplitPanel, Widget } from '@lumino/widgets';
import { ServiceManager, ServerConnection } from '@jupyterlab/services';
import { MathJaxTypesetter } from '@jupyterlab/mathjax2';
import { NotebookPanel, NotebookWidgetFactory, NotebookModelFactory } from '@jupyterlab/notebook';
import { CompleterModel, Completer, CompletionHandler, KernelConnector } from '@jupyterlab/completer';
import { editorServices } from '@jupyterlab/codemirror';
import { DocumentManager } from '@jupyterlab/docmanager';
import { DocumentRegistry } from '@jupyterlab/docregistry';
import { RenderMimeRegistry, standardRendererFactories as initialFactories } from '@jupyterlab/rendermime';
import { SetupCommands } from './commands';


//__webpack_public_path__ = "http://localhost:8082/packages/published/files/Notebooks/vnerozin/dist/";


//name: Main
export function main() {
    const settings = {
        baseUrl: 'http://localhost:8889',
        wsUrl: 'ws://localhost:8889',
        token: '77974bcd14909e7ee29330ac3657edc1e8d4ac39de425210',
        appUrl: '/',
        mathjaxUrl: 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js',
        mathjaxConfig: 'TeX-AMS_CHTML-full,Safe'
    };

    for (let key in settings)
          PageConfig.setOption(key, settings[key]);

    let _settings = ServerConnection.defaultSettings;

    _settings.baseUrl = PageConfig.getBaseUrl();
    _settings.appUrl = PageConfig.getOption('appUrl');
    _settings.wsUrl = PageConfig.getWsUrl();
    _settings.token = PageConfig.getToken();

    const manager = new ServiceManager({
        serverSettings: _settings
    });
    void manager.ready.then(() => {
        openNotebook(manager, 'test.ipynb');
    });
}


function openNotebook(manager, notebookPath) {
    // Initialize the command registry with the bindings.
    const commands = new CommandRegistry();
    const useCapture = true;

    // Setup the keydown listener for the document.
    document.addEventListener(
        'keydown',
        event => {
            commands.processKeydownEvent(event);
        },
        useCapture
    );

    const rendermime = new RenderMimeRegistry({
        initialFactories: initialFactories,
        latexTypesetter: new MathJaxTypesetter({
            url: PageConfig.getOption('mathjaxUrl'),
            config: PageConfig.getOption('mathjaxConfig')
        })
    });

    const opener = { open: (widget) => {} };
    const docRegistry = new DocumentRegistry();
    const docManager = new DocumentManager({
        registry: docRegistry,
        manager,
        opener
    });
    const mFactory = new NotebookModelFactory({});
    const editorFactory = editorServices.factoryService.newInlineEditor;
    const contentFactory = new NotebookPanel.ContentFactory({ editorFactory });

    const wFactory = new NotebookWidgetFactory({
        name: 'Notebook',
        modelName: 'notebook',
        fileTypes: ['notebook'],
        defaultFor: ['notebook'],
        preferKernel: true,
        canStartKernel: true,
        rendermime,
        contentFactory,
        mimeTypeService: editorServices.mimeTypeService
    });
    docRegistry.addModelFactory(mFactory);
    docRegistry.addWidgetFactory(wFactory);

    const nbWidget = docManager.open(notebookPath);
    const palette = new CommandPalette({ commands });
    palette.addClass('notebookCommandPalette');

    const editor = nbWidget.content.activeCell && nbWidget.content.activeCell.editor;
    const model = new CompleterModel();
    const completer = new Completer({ editor, model });
    const sessionContext = nbWidget.context.sessionContext;
    const connector = new KernelConnector({
        session: sessionContext.session
    });
    const handler = new CompletionHandler({ completer, connector });

    void sessionContext.ready.then(() => {
        handler.connector = new KernelConnector({
            session: sessionContext.session
        });
    });

    // Set the handler's editor.
    handler.editor = editor;

    // Listen for active cell changes.
    nbWidget.content.activeCellChanged.connect((sender, cell) => {
        handler.editor = cell && cell.editor;
    });

    // Hide the widget when it first loads.
    completer.hide();

    const panel = new SplitPanel();
    panel.id = 'main';
    panel.orientation = 'horizontal';
    panel.spacing = 0;
    SplitPanel.setStretch(palette, 0);
    SplitPanel.setStretch(nbWidget, 1);
    panel.addWidget(palette);
    panel.addWidget(nbWidget);

    let view = DG.View.create();
    view.name = 'Notebook';
    grok.shell.addView(view);

    Widget.attach(panel, view.root);
    Widget.attach(completer, view.root);

    // Handle resize events.
    view.root.addEventListener('resize', () => {
        panel.update();
    });

    view.root.classList.add('grok-notebook-view');

    SetupCommands(commands, palette, nbWidget, handler);

    console.log('Example started!');
}
