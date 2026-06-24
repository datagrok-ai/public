import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import '../css/notebooks.css';

import {PageConfig} from '@jupyterlab/coreutils';
import {CommandRegistry} from '@lumino/commands';
import {ServerConnection, ServiceManager} from '@jupyterlab/services';
import {MathJaxTypesetter} from '@jupyterlab/mathjax2';
import {
  CellTypeSwitcher,
  NotebookActions,
  NotebookModelFactory,
  NotebookPanel,
  NotebookWidgetFactory,
} from '@jupyterlab/notebook';
import {Completer, CompleterModel, CompletionHandler, KernelConnector} from '@jupyterlab/completer';
import {editorServices} from '@jupyterlab/codemirror';
import {DocumentManager} from '@jupyterlab/docmanager';
import {DocumentRegistry} from '@jupyterlab/docregistry';
import {RenderMimeRegistry, standardRendererFactories as initialFactories} from '@jupyterlab/rendermime';
import {SetupCommands} from './commands';
import {editNotebook, setupEnvironment, getAuthToken} from './utils';

export const _package = new DG.Package();

class NotebookView extends DG.ViewBase {
  constructor(params, path) {
    super(params, path);

    this.TYPE = 'Notebook';
    this.PATH = '/notebook';

    this.openAsHtmlIcon = ui.bigButton('HTML', () => {
      this.htmlMode().then();
    }, 'Open as HTML');
    this.openAsHtmlIcon.classList.add('notebooks-action-button');
    this.editIcon = ui.bigButton('EDIT', () => {
      this.editMode().then();
    }, 'Edit notebook');
    this.editIcon.classList.add('notebooks-action-button');

    this.html = '';
    this.saveAsComboPopup = ui.comboPopup(ui.iconFA('arrow-to-bottom', () => {
    }, 'Download'), ['As HTML', 'As PDF'], (item) => {
      if (item === 'As HTML') {
        const a = document.createElement('a');
        a.download = `${this.notebook.name}.html`;
        a.href = URL.createObjectURL(new Blob([this.html], {type: 'text/plain'})).toString();
        a.click();
      } else {
        const w = window.open('', '', 'height=650,width=900,top=100,left=150');
        const d = w['document'];
        d['body']['innerHTML'] = this.html;
        d['title'] = this.notebook.name;
        w.print();
      }
    });

    if (params !== null) {
      this.id = 'id' in params ? params['id'] :
        ((path !== null && path !== undefined) ? path.replace(`${this.PATH}/`, '') :
          null);

      const edit = 'edit' in params ? params['edit'] : false;
      const html = 'html' in params ? params['html'] : null;

      if (edit)
        this.editMode().then();
      else
        this.htmlMode(html).then();
    }
  }

  get type() {
    return this.TYPE;
  };

  get helpUrl() {
    return '/help/compute/jupyter-notebook.md';
  }

  get path() {
    return `${this.PATH}/${this.id}`;
  };

  get name() {
    return (this.notebook !== null && this.notebook !== undefined) ? this.notebook.name : 'Notebook';
  };

  set name(s) {
    super.name = s;
  }

  get entity() {
    return (this.notebook !== null && this.notebook !== undefined) ? this.notebook.d : null;
  }

  set entity(e) {
    super.entity = e;
  }

  getIcon() {
    const img = document.createElement('img');
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
    ui.empty(this.root);
    if (this.html === null)
      this.html = await toHtml(this.notebook);
    const iframe = document.createElement('iframe');
    const blob = new Blob([this.html.replace(/\u00a0/g, ' ')], {type: 'text/html'});
    iframe.src = URL.createObjectURL(blob);
    iframe.classList.add('grok-notebook-view-iframe');
    const sheet = new CSSStyleSheet();
    sheet.replaceSync('.grok-notebook-view-iframe {width: 100%; height:100%; border: none;display: flex;flex-grow: 1;}');

    const container = ui.div(null, 'grok-notebook-view-container');
    const shadowRoot = container.attachShadow({'mode': 'open'});
    shadowRoot.appendChild(iframe);
    shadowRoot.adoptedStyleSheets = [sheet];
    const view = ui.div([container], 'd4-root,d4-flex-col');
    this.setRibbonPanels([[this.saveAsComboPopup, this.editIcon]], true);
    this.editIcon.parentNode.parentNode.style.flexShrink = '0';
    this.root.appendChild(view);
  }

  async getEnvironmentsInput() {
    // TODO: Deprecate 'default' tricks after default environment implementation
    this.envs = [DG.ScriptEnvironment.create('default')];
    this.envs.push(...(await grok.dapi.environments.list()));
    let selected = this.envs.filter((e) => e.name === this.notebook.name);
    selected = (selected.length === 0) ? this.envs[0].name : selected[0].name;
    const environmentInput = ui.input.choice('Environment:', {value: selected, items: this.envs.map((e) => e.name)});
    environmentInput.onChanged.subscribe(async (value) => {
      const environment = this.envs.filter((e) => e.name === value)[0];
      if (this.notebook.environment === 'python3' && environment.name === 'default')
        return;
      if (this.notebook.environment !== environment.name) {
        this.notebook = await grok.dapi.notebooks.find(this.notebook.id);
        this.notebook.environment = environment.name;
        if (environment.name !== 'default') {
          ui.setUpdateIndicator(this.root, true);
          await setupEnvironment(environment, CONTAINER_ID);
          ui.setUpdateIndicator(this.root, false);
        }
        this.editMode();
      }
    });
    return environmentInput;
  }

  async editMode() {
    await this.initNotebook();
    ui.empty(this.root);

    this.subs.push(grok.events.onEvent('d4-entity-edited').subscribe((e) => {
      if (this.notebook !== undefined && this.notebook !== null && e.id === this.notebook.id)
        this.name = e.name;
    }));

    const notebookPath = await editNotebook(this.notebook, CONTAINER_ID);
    const manager = new ServiceManager({serverSettings: NotebookView.getSettings()});
    await manager.ready;

    const renderMime = new RenderMimeRegistry({
      initialFactories: initialFactories,
      latexTypesetter: new MathJaxTypesetter({
        url: PageConfig.getOption('mathjaxUrl'),
        config: PageConfig.getOption('mathjaxConfig'),
      }),
    });

    const opener = {
      open: (_) => {
      },
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
      mimeTypeService: editorServices.mimeTypeService,
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

    const iframe = document.createElement('iframe');
    iframe.classList.add('notebooks-edit-iframe');

    const commands = new CommandRegistry();
    const useCapture = true;

    // Datagrok detaches a view's root on tab switch; re-attaching the iframe forces the browser
    // to reload it to a blank document. The notebook model lives in memory, so we rebuild the
    // iframe document and re-attach the existing widget on every load to preserve all cell data.
    const mountNotebook = () => {
      const iframeDocument = iframe.contentDocument || iframe.contentWindow.document;
      iframeDocument.open();
      iframeDocument.write('<html><head></head><body></body></html>');
      iframeDocument.close();
      function addLink(path) {
        const link = document.createElement('link');
        link.rel = 'stylesheet';
        link.type = 'text/css';
        link.href = _package.webRoot + 'dist/' + path;
        iframeDocument.head.append(link);
      }
      addLink('styles/jupyter-styles.css');
      iframeDocument.documentElement.classList.add('notebooks-edit-html');
      iframeDocument.body.classList.add('notebooks-edit-body');
      const container = iframeDocument.createElement('div');
      container.classList.add('notebooks-edit-container');
      iframeDocument.body.append(container);
      iframeDocument.addEventListener('keydown', (event) => {
        commands.processKeydownEvent(event);
      }, useCapture);
      // iframe.addEventListener('mousemove', _event => {
      //   iframe.contentWindow.focus();
      // }, useCapture);
      // container.append(nbWidget.node);
      // container.append(completer.node);
      // MessageLoop.sendMessage(nbWidget, Widget.Msg.BeforeAttach);
      container.appendChild(nbWidget.node);
      // MessageLoop.sendMessage(nbWidget, Widget.Msg.AfterAttach);
      // MessageLoop.sendMessage(completer, Widget.Msg.BeforeAttach);
      container.appendChild(completer.node);
      // MessageLoop.sendMessage(completer, Widget.Msg.AfterAttach);
    };

    iframe.addEventListener('load', mountNotebook);
    this.root.appendChild(iframe);
    mountNotebook();

    handler.editor = editor;
    nbWidget.content.activeCellChanged.connect((sender, cell) => {
      handler.editor = cell !== null && cell !== undefined ? cell.editor : null;
      if (cell && cell.editor) cell.editor.focus();
    });

    nbWidget.content.node.addEventListener('click', (event) => {
      const cells = nbWidget.content.widgets; // all cells in the notebook
      for (let i = 0; i < cells.length; i++) {
        if (cells[i].node.contains(event.target)) {
          nbWidget.content.activeCellIndex = i; // sets the active cell by index
          break;
        }
      }
    });

    completer.hide();


    // if (this.environmentInput === undefined)
    //   this.environmentInput = await this.getEnvironmentsInput();

    this.setRibbonPanels([
      [
        this.openAsHtmlIcon,
        ui.iconFA('file-export', () => {
          grok.shell.addView(DG.ScriptView.create(DG.Script.create(this.notebookToCode(nbWidget.content))));
        }, 'Open as script'),
      ],
      [
        ui.iconFA('save', () => {
          nbWidget.context.save();
        }, 'Save notebook'),
        ui.iconFA('plus', () => NotebookActions.insertBelow(nbWidget.content), 'Insert a cell before'),
        ui.iconFA('cut', () => NotebookActions.cut(nbWidget.content), 'Cut cell'),
        ui.iconFA('copy', () => NotebookActions.copy(nbWidget.content), 'Copy cell'),
        ui.iconFA('paste', () => NotebookActions.copy(nbWidget.content), 'Paste cell'),
        ui.iconFA('play', () => NotebookActions.runAndAdvance(nbWidget.content, nbWidget.context.sessionContext), 'Run cell'),
        ui.iconFA('stop', () => nbWidget.context.sessionContext.session.kernel.interrupt(), 'Interrupt Kernel'),
        ui.iconFA('redo', () => {
          ui.dialog({title: 'Restart Kernel'})
            .onOK(() => sessionContext.restartKernel())
            .show();
        }, 'Restart'),
        ui.iconFA('forward', () => {
          ui.dialog({title: 'Restart Kernel and run all cells'})
            .onOK(() => sessionContext.restartKernel().then((restarted) => {
              if (restarted) NotebookActions.runAll(nbWidget.content, nbWidget.context.sessionContext);
            }))
            .show();
        }, 'Run all'),
        new CellTypeSwitcher(nbWidget.content).node,
      ],
      // [this.environmentInput.root],
    ], true);
    nbWidget.toolbar.hide();
    this.openAsHtmlIcon.parentNode.parentNode.style.flexShrink = '0';

    this.root.classList.add('grok-notebook-view');

    SetupCommands(commands, nbWidget, handler);
  }

  notebookToCode(jnb) {
    const inputRegex = /(.*) = api\.tables\.download/g;
    const outputRegex = /(.*) = api\.tables\.upload\((.*)\)/g;
    const script = [
      `#name: ${this.notebook.name}\n`,
      `#description: ${this.notebook.description}\n`,
      '#language: python\n',
      '#tags: notebook\n',
    ];
    const inputs = [];
    const outputs = [];
    const body = ['\n'];
    jnb.widgets.forEach((cell) => {
      const text = cell.model.value.text;
      if (cell.model.type === 'markdown' || cell.model.type === 'raw') {
        const lines = text.split('\n');
        for (let n = 0; n < lines.length; n++)
          lines[n] = `#${lines[n]}\n`;
        body.push(...lines);
        body.push('\n');
      } else if (cell.model.type === 'code') {
        const lines = text.split('\n').filter((l) => !l.startsWith('%'));
        for (const line of lines) {
          let match = inputRegex.exec(line);
          if (match !== null)
            inputs.push(`#input: dataframe ${match[1]}\n`);
          match = outputRegex.exec(line);
          if (match !== null)
            outputs.push(`#output: dataframe ${match[1]}\n`);
        }
        body.push(lines.join('\n'));
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
      baseUrl: 'http://localhost:8080/notebook', // host will be removed, needed only to extract path and use it in proxy request to container
      mathjaxUrl: 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js',
      mathjaxConfig: 'TeX-AMS_CHTML-full,Safe',
    };

    for (const key in _settings) PageConfig.setOption(key, _settings[key]);

    const settings = ServerConnection.defaultSettings;

    settings.fetch = async (info, _) => {
      const url = new URL(info.url);
      const path = url.pathname + url.search;
      const params = {method: info.method, headers: info.headers};
      if (info.body && info.body instanceof ReadableStream) {
        const reader = info.body.getReader();
        const chunks = [];
        let totalLength = 0;

        while (true) {
          const {done, value} = await reader.read();
          if (done) break;

          chunks.push(value);
          totalLength += value.length;
        }
        const bodyBuffer = new Uint8Array(totalLength);
        let offset = 0;
        for (const chunk of chunks) {
          bodyBuffer.set(chunk, offset);
          offset += chunk.length;
        }

        params.body = bodyBuffer.buffer;
      }
      return await grok.dapi.docker.dockerContainers.fetchProxy(CONTAINER_ID, path, params);
    };

    settings.WebSocket = DockerWebSocket;

    settings.baseUrl = PageConfig.getBaseUrl();
    settings.wsUrl = PageConfig.getWsUrl();
    settings.token = PageConfig.getToken();

    return settings;
  }
}

function createDockerWebSocket(url, _) {
  const uri = new URL(url);
  return grok.dapi.docker.dockerContainers.webSocketProxySync(CONTAINER_ID, uri.pathname + uri.search);
}

const DockerWebSocket = new Proxy(createDockerWebSocket, {
  construct(target, args) {
    const ws = target(...args);
    let connected = false;
    const buffer = [];
    return new Proxy(ws, {
      set(target, prop, value) {
        if (prop === 'onmessage') {
          target[prop] = function(event) {
            if (event.data === 'CONNECTED') {
              connected = true;
              if (buffer.length !== 0) {
                for (const m of buffer)
                  ws.send(m);
                buffer.length = 0;
              }
              return;
            }
            value(event);
          };
        } else
          target[prop] = value;

        return true;
      },

      get(target, prop) {
        if (prop === 'send') {
          return new Proxy(target[prop], {
            apply: (sendMethod, thisArg, argumentsList) => {
              if (typeof argumentsList[0] === 'string')
                argumentsList[0] = argumentsList[0].replace(/@USER_API_KEY/g, `'${SESSION_TOKEN}'`);
              if (!connected)
                buffer.push(argumentsList[0]);
              else
                return Reflect.apply(sendMethod, target, argumentsList);
            },
          });
        } else if (prop === 'close') {
          // Intercept close method to ensure correct context
          return function(...args) {
            console.log('Calling WebSocket close...');
            return Reflect.apply(target[prop], target, args); // Call close on the actual WebSocket
          };
        } else
          return Reflect.get(target, prop);
      },
    });
  },
});

let CONTAINER_ID;
let SESSION_TOKEN;

async function toHtml(notebook) {
  const arrayBuffer = await convertNotebook(JSON.stringify(notebook.notebook));
  return new TextDecoder('utf-8').decode(arrayBuffer);
}

//tags: init
export async function initContainer() {
  const container = await grok.dapi.docker.dockerContainers.filter('Notebooks-jupyter-notebook').first();
  CONTAINER_ID = container.id;
  if (container.status !== 'started' && !container.status.startsWith('pending') && container.status !== 'checking') {
    const progress = DG.TaskBarProgressIndicator.create('Starting Jupyter Notebook...');
    try {
      await grok.dapi.docker.dockerContainers.run(CONTAINER_ID, true);
      // wait additional time, because nginx tends to start faster
      await new Promise((res, _) => setTimeout(() => res(), 3000));
    } finally {
      progress.close();
    }
  }
  SESSION_TOKEN = getAuthToken();
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

//name: convertNotebook
//description: Converts notebook file content into specified format
//input: string notebook
//input: string format
//input: bool execute
//output: blob result
export async function convertNotebook(notebook, format = 'html', execute = false) {
  notebook = notebook.replace(/@USER_API_KEY/g, `'${SESSION_TOKEN}'`);
  const response = await grok.dapi.docker.dockerContainers.fetchProxy(CONTAINER_ID, `/notebook/helper/notebooks/convert?format=${format}&execute=${execute}`,
    {method: 'POST', body: new TextEncoder().encode(notebook)});
  if (response.status > 201)
    throw response.statusText;
  return response.arrayBuffer();
}
