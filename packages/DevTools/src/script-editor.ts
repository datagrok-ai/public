import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {applyCodeMirror} from "./utils/code-mirror-check";

export function initScriptEditor(view: DG.View): void {
  applyCodeMirror(view, (editor) => {
    const doc = editor.getDoc();
    const lineCount = doc.lineCount();
    for (let i = 0; i < lineCount; i++) {
      if (doc.getLine(i).trim() == '//language: javascript') {
        setScriptRibbon(view, doc);
        editor.setCursor({line: 6, ch: 0});
        break;
        //editor.setFocus();
      }
    }
  });
}

function setScriptRibbon(v: DG.View, doc: any) {
  const viewsMenu = () => {
    const addView = (item) => {
      const cursor = doc.getCursor();
      const template = `//name: ${item}\n//language: javascript\nlet windows = grok.shell.windows;\nwindows.showToolbox = false;\nwindows.showProperties = false;\nwindows.showHelp = false;\n`;
      switch (item) {
      case 'Simple View':
        doc.replaceRange('grok.shell.newView(\'New view\',[])\n', cursor);
        break;
      case 'Table View':
        doc.replaceRange('grok.shell.addTableView(table)\n', cursor);
        break;
      case 'Viewer':
        doc.replaceRange('DG.Viewer.fromType(\'Scatter Plot\', table)\n', cursor);
        break;
      case 'Fixed':
        doc.replaceRange('grok.shell.dockManager.dock(ui.divText(\'Dock Content\'), \'right\', null, \'Title\');\n', cursor);
        break;
      case 'Floating':
        doc.replaceRange('grok.shell.dockManager.dock(ui.divText(\'Dock Content\')).container.float();\n', cursor);
        break;
      case 'Gallery Grid':
        doc.setValue(template + `\nlet view = grok.shell.newView('Gallery grid');view.box = true;\nlet gallery = ui.divH([],'grok-gallery-grid');\nview.append(gallery);`);
        break;
      case 'Flex Box':
        doc.setValue(template + `\nlet view = grok.shell.newView('Flex box' [\n  ui.divV([\n    ui.divH([]),\n  ])\n]);`);
        break;
      case 'Tab Navbar':
        doc.setValue(template + `\nlet view = grok.shell.newView('Tab navbar' [\n  ui.tabControl({\n    'First': ui.panel('First tab'),\n    'Second': ui.panel('Second tab'),\n  })\n]);`);
        break;
      case '2-Column':
        doc.setValue(template + `\nlet view = grok.shell.newView('Splitter 2-Column' [\n  ui.splitH([\n    ui.panel('Column 1'),\n    ui.panel('Column 2'),\n  ])\n]);\nview.box = true;`);
        break;
      case '2-Row':
        doc.setValue(template + `\nlet view = grok.shell.newView('Splitter 2-Row' [\n  ui.splitV([\n    ui.panel('Row 1'),\n    ui.panel('Row 2'),\n  ])\n]);\nview.box = true;`);
        break;
      case '1-Column, 2-Row':
        doc.setValue(template + `\nlet view = grok.shell.newView('Splitter 1-Column, 2-Row' [\n  ui.splitH([\n    ui.panel('Collumn 1'),\n    ui.splitV([\n      ui.panel('Row 1'),\n      ui.panel('Row 2'),\n    ]),\n  ])\n]);\nview.box = true;`);
        break;
      case '2-Column, 2-Row':
        doc.setValue(template + `\nlet view = grok.shell.newView('Splitter 2-Column, 2-Row' [\n  ui.splitV([\n    ui.splitH([\n      ui.panel('Row 1 Column 1'),\n      ui.panel('Row 1 Column 2'),\n    ]),\n    ui.splitH([\n      ui.panel('Row 2 Column 1'),\n      ui.panel('Row 2 Column 2'),\n    ]),\n]);\nview.box = true;`);
        break;
      }
    };

    const viewsPopup = DG.Menu.popup()
      .items(['Simple View', 'Table View'], addView)
      .separator()
      .items(['Viewer'], addView)
      .group('Dock Manager')
      .items(['Fixed', 'Floating'], addView)
      .endGroup()
      .separator()
      .group('Templates')
      .items(['Gallery Grid', 'Flex Box', 'Tab Navbar'], addView)
      .separator()
      .items(['2-Column', '2-Row', '1-Column, 2-Row', '2-Column, 2-Row'], addView)
      .endGroup();

    try {
      // @ts-ignore
      viewsPopup.find('Simple View')._check.appendChild(ui.iconFA('window-maximize'));
      // @ts-ignore
      viewsPopup.find('Table View')._check.appendChild(ui.iconFA('table'));
      // @ts-ignore
      viewsPopup.find('Viewer')._check.appendChild(ui.iconFA('chart-bar'));
      // @ts-ignore
      viewsPopup.group('Dock Manager').find('Fixed').parentMenu._check.appendChild(ui.iconFA('window-restore'));
    } catch (x) {

    } finally {
      viewsPopup.show();
    }
  };
  const layoutsMenu = () => {
    const addLayout = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'Panel':
        doc.replaceRange('ui.panel([ ])\n', cursor);
        break;
      case 'Div':
        doc.replaceRange('ui.div([ ])\n', cursor);
        break;
      case 'Box':
        doc.replaceRange('ui.box()\n', cursor);
        break;
      case 'Div-H':
        doc.replaceRange('ui.divH([ ])\n', cursor);
        break;
      case 'Div-V':
        doc.replaceRange('ui.divV([ ])\n', cursor);
        break;
      case '100%':
        doc.replaceRange('ui.block([ ])\n', cursor);
        break;
      case '25%':
        doc.replaceRange('ui.block25([ ])\n', cursor);
        break;
      case '50%':
        doc.replaceRange('ui.block50([ ])\n', cursor);
        break;
      case '75%':
        doc.replaceRange('ui.block75([ ])\n', cursor);
        break;
      case 'Split-H':
        doc.replaceRange('ui.splitH([ ])\n', cursor);
        break;
      case 'Split-V':
        doc.replaceRange('ui.splitV([ ])\n', cursor);
        break;
      case 'Inputs':
        doc.replaceRange('ui.inputs([ ])\n', cursor);
        break;
      case 'Button inputs':
        doc.replaceRange('ui.buttonsInput([ ])\n', cursor);
        break;
      }
    };

    const layoutsPopup = DG.Menu.popup()
      .group('Containers')
      .items(['Panel', 'Div', 'Box'], addLayout)
      .endGroup()
      .group('Flex box')
      .items(['Div-H', 'Div-V'], addLayout)
      .endGroup()
      .group('Blocks')
      .items(['100%', '75%', '50%', '25%'], addLayout)
      .endGroup()
      .group('Splitters')
      .items(['Split-H', 'Split-V'], addLayout)
      .endGroup()
      .group('Form')
      .items(['Inputs', 'Button inputs'], addLayout)
      .endGroup();

    layoutsPopup.show();
  };
  const elementsMenu = () => {
    const addElementText = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'H1':
        doc.replaceRange('ui.h1(\'Headline H1\')\n', cursor);
        break;
      case 'H2':
        doc.replaceRange('ui.h2(\'Headline H2\')\n', cursor);
        break;
      case 'H3':
        doc.replaceRange('ui.h3(\'Headline H3\')\n', cursor);
        break;
      case 'Paragraph':
        doc.replaceRange('ui.p(\'Paragraph\')\n', cursor);
        break;
      case 'Span':
        doc.replaceRange('ui.span([])\n', cursor);
        break;
      case 'Label':
        doc.replaceRange('ui.label(\'Label text\')\n', cursor);
        break;
      case 'Link':
        doc.replaceRange('ui.link(\'Text\',()=>{},\'Tooltip message\')\n', cursor);
        break;
      case 'Div Text':
        doc.replaceRange('ui.divText(\'Text block\')\n', cursor);
        break;
      case 'Inline Text':
        doc.replaceRange('ui.inlineText([\'Text line 1\',\'Text line 2\'])\n', cursor);
        break;
      }
    };
    const addElementButton = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'Button':
        doc.replaceRange('ui.button(\'name\',()=>{})\n', cursor);
        break;
      case 'Big Button':
        doc.replaceRange('ui.bigButton(\'name\',()=>{})\n', cursor);
        break;
      case 'Icon Button':
        doc.replaceRange('ui.button(ui.iconFA(\'smile\'),()=>{})\n', cursor);
        break;
      }
    };
    const addElementInput = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'Int':
        doc.replaceRange('ui.input.int(\'Label\', {value: 0})\n', cursor);
        break;
      case 'String':
        doc.replaceRange('ui.input.string(\'Label\', {value: \'Value\'})\n', cursor);
        break;
      case 'Date':
        doc.replaceRange('ui.input.date(\'Label\', {value: dayjs(\'1970-05-10\')})\n', cursor);
        break;
      case 'Bool':
        doc.replaceRange('ui.input.bool(\'Label\', {value: true})\n', cursor);
        break;
      case 'Choice':
        doc.replaceRange('ui.input.choice(\'Label\', {value: \'Val 1\', items: [\'Val 1\', \'Val 2\', \'Val 3\']})\n', cursor);
        break;
      case 'Multi Choice':
        doc.replaceRange('ui.input.multiChoice(\'Label\', {value: [\'Val 1\', \'Val 2\'], items: [\'Val 1\', \'Val 2\', \'Val 3\']})\n', cursor);
        break;
      case 'Text Area':
        doc.replaceRange('ui.input.textArea(\'Label\', {value: \'Value\'})\n', cursor);
        break;
      case 'Switch':
        doc.replaceRange('ui.input.toggle(\'demo\', {value: false})\n', cursor);
        break;
      case 'Column':
        doc.replaceRange('ui.input.column(\'Label\', {table: table, value: table.col(\'age\')})\n', cursor);
        break;
      case 'Columns':
        doc.replaceRange('ui.input.columns(\'Label\', {table: table})\n', cursor);
        break;
      case 'Table':
        doc.replaceRange('ui.input.table(\'Label\', {value: tables[0], items: tables, onValueChanged: (value) => grok.shell.info(value.name))\n', cursor);
        break;
      case 'Molecule':
        doc.replaceRange('ui.input.molecule(\'Label\', {value: \'CN1CCC(O)(CC1)c2ccccc2\'})\n', cursor);
        break;
      }
    };
    const addElementContent = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'HTML Table':
        doc.replaceRange('ui.tableFromMap({\'Key 1\':\'Val 1\', \'Key 2\':\'Val 2\'})\n', cursor);
        break;
      case 'List':
        doc.replaceRange('ui.list([\'Line 1\',\'Line 2\'])\n', cursor);
        break;
      case 'Image':
        doc.replaceRange('ui.image(\'src...\',100,100,{target:\'url\'})\n', cursor);
        break;
      case 'IconFA':
        doc.replaceRange('ui.iconFA(\'search\',()=>{})\n', cursor);
        break;
      case 'Add':
        doc.replaceRange('ui.icons.add(() => {}, \'Add\')\n', cursor);
        break;
      case 'Close':
        doc.replaceRange('ui.icons.close(() => {}, \'Close\')\n', cursor);
        break;
      case 'Copy':
        doc.replaceRange('ui.icons.copy(() => {}, \'Copy\')\n', cursor);
        break;
      case 'Delete':
        doc.replaceRange('ui.icons.delete(() => {}, \'Delete\')\n', cursor);
        break;
      case 'Edit':
        doc.replaceRange('ui.icons.edit(() => {}, \'Edit\')\n', cursor);
        break;
      case 'Filter':
        doc.replaceRange('ui.icons.filter(() => {}, \'Filter\')\n', cursor);
        break;
      case 'Help':
        doc.replaceRange('ui.icons.help(() => {}, \'Help\')\n', cursor);
        break;
      case 'Info':
        doc.replaceRange('ui.icons.info(() => {}, \'Info\')\n', cursor);
        break;
      case 'Remove':
        doc.replaceRange('ui.icons.remove(() => {}, \'Remove\')\n', cursor);
        break;
      case 'Save':
        doc.replaceRange('ui.icons.save(() => {}, \'Save\')\n', cursor);
        break;
      case 'Search':
        doc.replaceRange('ui.icons.search(() => {}, \'Search\')\n', cursor);
        break;
      case 'Settings':
        doc.replaceRange('ui.icons.settings(() => {}, \'Settings\')\n', cursor);
        break;
      case 'Sync':
        doc.replaceRange('ui.icons.sync(() => {}, \'Sync\')\n', cursor);
        break;
      case 'Undo':
        doc.replaceRange('ui.icons.undo(() => {}, \'Undo\')\n', cursor);
        break;
      }
    };

    const elementsPopup = DG.Menu.popup()
      .group('Typography')
      .items(['H1', 'H2', 'H3', 'Paragraph', 'Span', 'Label', 'Link', 'Div Text', 'Inline Text'], addElementText)
      .endGroup()
      .group('Inputs')
      .items(['Int', 'String', 'Date', 'Bool'], addElementInput)
      .separator()
      .items(['Choice', 'Multi Choice', 'Text Area', 'Switch'], addElementInput)
      .separator()
      .items(['Column', 'Columns', 'Table', 'Molecule'], addElementInput)
      .endGroup()
      .group('Buttons')
      .items(['Button', 'Big Button', 'Icon Button'], addElementButton)
      .endGroup()
      .group('Content')
      .items(['HTML Table', 'List', 'Image', 'IconFA'], addElementContent)
      .endGroup()
      .group('Special Icons')
      .items(['Add', 'Close', 'Copy', 'Delete', 'Edit', 'Filter', 'Help', 'Info', 'Remove', 'Save', 'Search', 'Settings', 'Sync', 'Undo'], addElementContent)
      .endGroup();

    try {
      const add = (name: string, makeIcon: Function) => {
        elementsPopup.group('Special Icons').find(name)._check.append(makeIcon(() => { }, ''));
      }

      add('Add', ui.icons.add);
      add('Close', ui.icons.close);
      add('Copy', ui.icons.copy);
      add('Delete', ui.icons.delete);
      add('Edit', ui.icons.edit);
      add('Filter', ui.icons.filter);
      add('Help', ui.icons.help);
      add('Info', ui.icons.info);
      add('Remove', ui.icons.remove);
      add('Save', ui.icons.save);
      add('Search', ui.icons.search);
      add('Settings', ui.icons.settings);
      add('Sync', ui.icons.sync);
      add('Undo', ui.icons.undo);
    } catch (x) {
      console.log(x);
    } finally {
      elementsPopup.show();
    }
  };
  const componentsMenu = () => {
    const addComponent = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'Accordion':
        doc.replaceRange('let acc = ui.accordion();\nacc.addPane(\'Accordion Pane\', () => ui.divV([\n  ui.p(\'Pane content\')\n  ]));', cursor);
        break;
      case 'Combo-popup':
        doc.replaceRange('ui.comboPopup(\'Label\', [\'Value 1\', \'Value 2\'], (item) => {})\n', cursor);
        break;
      case 'Column Combo-popup':
        doc.replaceRange('DG.ColumnComboBox.create(grok.data.demo.demog(), (c) => c.type === \'string\')\n', cursor);
        break;
      case 'Iframe':
        doc.replaceRange('ui.iframe( { src: \'https://en.m.wikipedia.org/wiki\', width: \'400\', height: \'200\' })\n', cursor);
        break;
      case 'Info Bar':
        doc.replaceRange('ui.info(\n  ui.divText(\'Content \'),\n  ui.h1(\'Title optional\'),\n  false\n),', cursor);
        break;
      case 'Markdown':
        doc.replaceRange('ui.markdown(`# Title \n Content`)\n', cursor);
        break;
      case 'Menu':
        doc.replaceRange('grok.shell.topMenu\n.group(\'Manu name\')\n.item(\'Item\', ()=>{})\n.group(\'Group\')\n.item(\'Group item\',()=>{})\n.endGroup()\n.endGroup()', cursor);
        break;
      case 'Popup Menu':
        doc.replaceRange('let showMenu = () => {\n  DG.Menu.popup()\n    .item(\'Item 1\', ()=>{})\n    .show()\n};', cursor);
        break;
      case 'Sidebar':
        doc.replaceRange('grok.shell.sidebar.addPane(\'Name\', () => ui.divText(\'Panel content\'), ui.iconFA(\'universal-access\'));', cursor);
        break;
      case 'Tab Control':
        doc.replaceRange('ui.tabControl({\n  \'FIRST\': ui.panel([]),\n  \'SECOND\': ui.panel([])\n});', cursor);
        break;
      case 'Tag Editor':
        doc.replaceRange('let editor = DG.TagEditor.create();\neditor.addTag(\'Tag name\')', cursor);
        break;
      case 'Taskbar Progress':
        doc.replaceRange('let p = DG.TaskBarProgressIndicator.create(\'Progress...\');\nsetTimeout(() => {\n  p.close();\n}, 3000);', cursor);
        break;
      case 'Tooltip':
        doc.replaceRange('ui.tooltip.bind(ui.span([\'Element\']), \'Tooltip message\')\n', cursor);
        break;
      case 'Tree view':
        doc.replaceRange('let tree = ui.tree();\nlet group = tree.group(\'group 1\', 1);\ngroup.enableCheckBox();\ngroup.item(\'item 1.1\').enableCheckBox();\ngroup.item(\'item 1.2\');', cursor);
        break;
      }
    };

    const addComponentDialog = (item) => {
      const cursor = doc.getCursor();
      switch (item) {
      case 'Standard':
        doc.replaceRange('ui.dialog(\'Title\',\'help url ...\')\n  .add(ui.p(\'Dialog content\'))\n  .onOK(()=>{})\n  .show();', cursor);
        break;
      case 'Modal':
        doc.replaceRange('ui.dialog(\'Title\',\'help url ...\')\n  .add(ui.p(\'Dialog content\'))\n  .onOK(()=>{})\n  .showModal();', cursor);
        break;
      case 'Full Screen':
        doc.replaceRange('ui.dialog(\'Title\',\'help url ...\')\n  .add(ui.p(\'Dialog content\'))\n  .onOK(()=>{})\n  .showModal(true);', cursor);
        break;
      }
    };

    const componentsPopup = DG.Menu.popup()
      .items(['Accordion', 'Combo-popup', 'Column Combo-popup'], addComponent)
      .group('Dialogs')
      .items(['Standard', 'Modal', 'Full Screen'], addComponentDialog)
      .endGroup()
      .items(['Iframe', 'Info Bar', 'Markdown', 'Menu', 'Popup Menu', 'Sidebar', 'Tab Control', 'Tag Editor', 'Taskbar Progress', 'Tooltip', 'Tree view'], addComponent);

    componentsPopup.show();
  };

  const view = ui.iconFA('window-maximize', () => viewsMenu);
  view.className = 'grok-icon far fa-window-maximize';
  const viewBtn = ui.div([view], 'd4-combo-popup');
  viewBtn.addEventListener('click', viewsMenu);

  const layouts = ui.iconFA('columns', () => layoutsMenu);
  layouts.className = 'grok-icon far fa-columns';
  const layoutsBtn = ui.div([layouts], 'd4-combo-popup');
  layoutsBtn.addEventListener('click', layoutsMenu);

  const elements = ui.iconFA('sliders-h');
  elements.className = 'grok-icon far fa-sliders-h';
  const elementsBtn = ui.div([elements], 'd4-combo-popup');
  elementsBtn.addEventListener('click', elementsMenu);

  const components = ui.iconFA('cubes');
  components.className = 'grok-icon fas fa-cubes';
  const componentsBtn = ui.div([components], 'd4-combo-popup');
  componentsBtn.addEventListener('click', componentsMenu);

  const panels = v.getRibbonPanels();
  const newPanels = [
    ui.tooltip.bind(viewBtn, 'Views'),
    ui.tooltip.bind(layoutsBtn, 'Layouts'),
    ui.tooltip.bind(elementsBtn, 'Elements'),
    ui.tooltip.bind(componentsBtn, 'Components'),
  ];

  if (JSON.stringify(newPanels) !== JSON.stringify(panels[panels.length - 2])) {
    v.setRibbonPanels([
      ...panels,
      newPanels,
    ]);
  }
}
