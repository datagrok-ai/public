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

import {PageConfig} from '@jupyterlab/coreutils';
import {CommandRegistry} from '@lumino/commands';
import {Widget} from '@lumino/widgets';
import {ServiceManager, ServerConnection} from '@jupyterlab/services';
import {MathJaxTypesetter} from '@jupyterlab/mathjax2';
import {
  NotebookPanel,
  NotebookWidgetFactory,
  NotebookModelFactory,
  NotebookActions,
  CellTypeSwitcher
} from '@jupyterlab/notebook';
import {CompleterModel, Completer, CompletionHandler, KernelConnector} from '@jupyterlab/completer';
import {editorServices} from '@jupyterlab/codemirror';
import {DocumentManager} from '@jupyterlab/docmanager';
import {DocumentRegistry} from '@jupyterlab/docregistry';
import {RenderMimeRegistry, standardRendererFactories as initialFactories} from '@jupyterlab/rendermime';
import {SetupCommands} from './commands';
import {sessionContextDialogs} from '@jupyterlab/apputils';


class NotebookView extends DG.ViewBase {
  constructor(params, path) {
    super(params, path);

    this.TYPE = 'Notebook';
    this.PATH = '/notebook';

    this.openAsHtmlIcon = ui.bigButton('HTML', () => {
      this.htmlMode().then()
    }, 'Open as HTML');
    this.openAsHtmlIcon.style.backgroundColor = 'var(--green-2)';
    this.editIcon = ui.bigButton('EDIT', () => {
      this.editMode().then()
    }, 'Edit notebook');
    this.editIcon.style.backgroundColor = 'var(--green-2)';

    this.html = '';
    this.saveAsComboPopup = ui.comboPopup(ui.iconFA('arrow-to-bottom', () => {
    }), ['As HTML', 'As PDF'], (item) => {
      if (item === 'As HTML') {
        let a = document.createElement('a');
        a.download = `${this.notebook.name}.html`;
        a.href = URL.createObjectURL(new Blob([this.html], {type: 'text/plain'})).toString();
        a.click();
      } else {
        let w = window.open('', '', 'height=650,width=900,top=100,left=150');
        let d = w['document'];
        d['body']['innerHTML'] = this.html;
        d['title'] = this.notebook.name;
        w.print();
      }
    });

    if (params !== null) {
      this.id = 'id' in params ? params['id']
        : ((path !== null && path !== undefined) ? path.replace(`${this.PATH}/`, '')
          : null);

      let edit = 'edit' in params ? params['edit'] : false;
      let html = 'html' in params ? params['html'] : null;

      if (edit) this.editMode().then();
      else this.htmlMode(html).then();
    }
  }

  get type() {
    return this.TYPE
  };

  get helpUrl() {
    return '/help/compute/jupyter-notebook.md';
  }

  get path() {
    return `${this.PATH}/${this.id}`
  };

  get name() {
    return (this.notebook !== null && this.notebook !== undefined) ? this.notebook.name : 'Notebook'
  };

  set name(s) {
    super.name = s;
  }

  get entity() {
    return (this.notebook !== null && this.notebook !== undefined) ? this.notebook.d : null;
  }

  set entity(e) {
    super.entity = e
  }

  getIcon() {
    let img = document.createElement('img');
    img.src = '/images/entities/jupyter.png';
    img.height = 18;
    img.width = 18;
    return img;
  };

  saveStateMap() {
    return {'id': this.id};
  }

  loadStateMap(stateMap) {
    this.open(stateMap['id']);
  }

  handlePath(path) {
    this.open(path.replace(`${this.PATH}/`, ''));
  }

  acceptsPath(path) {
    return path.startsWith(this.PATH);
  }

  open(id) {
    this.id = id;
    this.htmlMode(null).then();
  }

  async initNotebook() {
    this.notebook = await grok.dapi.notebooks.find(this.id);
    this.name = this.notebook.name;
  }

  async htmlMode(html = null) {
    this.html = html;
    await this.initNotebook();
    removeChildren(this.root);
    if (this.html === null)
      this.html = await this.notebook.toHtml();
    let iframe = document.createElement('iframe');
    iframe.src = 'data:text/html;base64,' + btoa(this.html.replace(/\u00a0/g, " "));
    iframe.classList.add('grok-notebook-view-iframe');
    let container = ui.div([iframe], 'grok-notebook-view-container');
    let view = ui.div([container], 'd4-root,d4-flex-col');
    this.setRibbonPanels([[this.saveAsComboPopup, this.editIcon]], true);
    this.editIcon.parentNode.parentNode.style.flexShrink = '0';
    this.root.appendChild(view);
  }

  async getEnvironmentsInput() {
    // TODO: Deprecate 'default' tricks after default environment implementation
    this.envs = [DG.ScriptEnvironment.create('default')];
    this.envs.push(...(await grok.dapi.environments.list()));
    let selected = this.envs.filter(e => e.name === this.notebook.name);
    selected = (selected.length === 0) ? this.envs[0].name : selected[0].name;
    let environmentInput = ui.input.choice('Environment:', {value: selected, items: this.envs.map(e => e.name)});
    environmentInput.onChanged.subscribe(async (value) => {
      let environment = this.envs.filter(e => e.name === value)[0];
      if (this.notebook.environment === 'python3' && environment.name === 'default')
        return;
      if (this.notebook.environment !== environment.name) {
        this.notebook = await grok.dapi.notebooks.find(this.notebook.id);
        this.notebook.environment = environment.name;
        if (environment.name !== 'default') {
          ui.setUpdateIndicator(this.root, true);
          await environment.setup();
          ui.setUpdateIndicator(this.root, false);
        }
        this.editorMode();
      }
    });
    return environmentInput;
  }

  async editMode() {
    await this.initNotebook();
    removeChildren(this.root);

    this.subs.push(grok.events.onEvent('d4-entity-edited').subscribe((e) => {
      if (this.notebook !== undefined && this.notebook !== null && e.id === this.notebook.id)
        this.name = e.name;
    }));

    let notebookPath = await this.notebook.edit();
    const manager = new ServiceManager({serverSettings: NotebookView.getSettings()});
    await manager.ready;

    // Initialize the command registry with the bindings.
    // Setup the keydown listener for the document.
    const commands = new CommandRegistry();
    const useCapture = true;
    document.addEventListener('keydown', event => {
      commands.processKeydownEvent(event);
    }, useCapture);

    const renderMime = new RenderMimeRegistry({
      initialFactories: initialFactories,
      latexTypesetter: new MathJaxTypesetter({
        url: PageConfig.getOption('mathjaxUrl'),
        config: PageConfig.getOption('mathjaxConfig')
      })
    });

    const opener = {
      open: (widget) => {
      }
    };
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
      NotebookActions.runAll(nbWidget.content, nbWidget.context.sessionContext).then();
    });

    handler.editor = editor;

    nbWidget.content.activeCellChanged.connect((sender, cell) => {
      handler.editor = cell !== null && cell !== undefined ? cell.editor : null;
    });

    completer.hide();

    Widget.attach(nbWidget, this.root);
    Widget.attach(completer, this.root);

    if (this.environmentInput === undefined)
      this.environmentInput = await this.getEnvironmentsInput();

    this.setRibbonPanels([
      [
        this.openAsHtmlIcon,
        ui.iconFA('file-export', () => {
          grok.shell.addView(DG.ScriptView.create(DG.Script.create(this.notebookToCode(nbWidget.content))));
        }, 'Open as script')
      ],
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
            .then(restarted => {
              if (restarted) NotebookActions.runAll(nbWidget.content, nbWidget.context.sessionContext);
            }),
          'Restart Kernel and run all cells'),
        new CellTypeSwitcher(nbWidget.content).node,
      ],
      [this.environmentInput.root],
    ], true);
    nbWidget.toolbar.hide();
    this.openAsHtmlIcon.parentNode.parentNode.style.flexShrink = '0';

    this.root.classList.add('grok-notebook-view');

    SetupCommands(commands, nbWidget, handler);
  }

  notebookToCode(jnb) {
    let inputRegex = /(.*) = grok_read/g;
    let outputRegex = /grok\((.*)\)/g;
    let script = [
      `#name: ${this.notebook.name}\n`,
      `#description: ${this.notebook.description}\n`,
      '#language: python\n',
      '#tags: notebook\n'
    ];
    let inputs = [];
    let outputs = [];
    let body = ['\n'];
    jnb.widgets.forEach((cell) => {
      let text = cell.model.value.text;
      if (cell.model.type === 'markdown' || cell.model.type === 'raw') {
        let lines = text.split('\n');
        for (let n = 0; n < lines.length; n++)
          lines[n] = `#${lines[n]}\n`;
        body.push(...lines);
        body.push('\n');
      } else if (cell.model.type === 'code') {
        let lines = text.split('\n').filter(l => !l.startsWith('%'));
        for (let line of lines) {
          let match = inputRegex.exec(line);
          if (match !== null)
            inputs.push(`#input: dataframe ${match[1]}\n`);
          match = outputRegex.exec(line);
          if (match !== null)
            outputs.push(`#output: dataframe ${match[1]}\n`);
        }
        lines = lines.filter(l => !l.includes('grok')).join('\n');
        body.push(lines);
        body.push('\n\n');
      }
    });
    if (inputs.length > 0) script.push(...inputs);
    if (outputs.length > 0) script.push(...outputs);
    script.push(...body);
    return script.join('');
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
//input: map params
//input: string path
//tags: view
//output: view result
export function notebookView(params = null, path = '') {
  return new NotebookView(params, path);
}


function removeChildren(node) {
  while (node.firstChild)
    node.removeChild(node.firstChild);
}
