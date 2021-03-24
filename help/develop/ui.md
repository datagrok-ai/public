<!-- TITLE: Grok UI -->
<!-- SUBTITLE: -->

<!-- This is a user-centric view on the Datagrok applications UI development-->

# Layouts

 ## Starting point
 To start building a layout you need either a [View](#simple-view) (in most cases) or a [Dialog](#dialogs).
 Every View or Dialog is a [panel container](#panels)

```javascript
let v = grok.shell.newView('Demo View');
v.append(ui.h1('Hello World'));
let d = ui.dialog('Demo Dialog');
d.add(ui.h1('Hello World'));
d.show();
```

 ## Containers
This is a simple container. It can contain any elements, such as inputs, images, text, etc.
It doesn't have own height. The container height depends on its children.

![Container preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/container.jpg)

```javascript
ui.div([ui.h1('Header'), ui.p('Paragraph text'), 'just text', DG.Viewer.scatterPlot(grok.data.demo.demog())])
```
If you place a container within a [box container](#boxes) , it will inherit the box size and show scroll bars automatically.
 ## Boxes
This is a fixed-size container. It doesn't depend on children element sizes, but shrinks them to certain size.

![Box container preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/box-container.jpg)

```javascript
let d = ui.div();
for (let i = 0; i < 100; i++)
    d.append(ui.p('More Text'));
var box = ui.box(d);
$(box).css('border', 'solid');
ui.div([ui.h1('Header'), box])
```
 ## Panels
 The panel is a simple container similar to the [Containers](#Containers). It has full available wide and its height depends on its children. Also, panels have 10px paddings for all sides.

 ![Panel preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/panel.jpg)

 ```javascript
 ui.panel([ui.h1('Header'), ui.p('Paragraph text'), 'just text', DG.Viewer.scatterPlot(grok.data.demo.demog())])
 ```

 ## Blocks
 The block layout is divided into horizontal sections, which take on the full width of the available screen. Their screen height is determined by their inner content. The width of sections can also be set to the following predefined ratios:
 - 1 block: 100%
 - 2 blocks:
   - 50% and 50%
   - 75% and 25%
   - 25% and 75%
 - 3 blocks:
   - 2 x 25% and 50%
   - 50% and 2 x 25%
 - 4 blocks: 4 x 25%

 Use the block layout if you want to display section-based content by placing elements next to each other.

![Blocks preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/blocks.jpg)

 ```javascript
 ui.block([ui.h1('100% block width')]);
 //two blocks next to each other
 ui.block75([ui.h1('75% block width')]);  
 ui.block25([ui.h1('25% block width')]);
 //three blocks next ot each other
 ui.block50([ui.h1('50% block width')]);
 ui.block25([ui.h1('25% block width')]);
 ui.block25([ui.h1('25% block width')]);
 //four blocks next ot each other
 ui.block25([ui.h1('25% block width')]);
 ui.block25([ui.h1('25% block width')]);
 ui.block25([ui.h1('25% block width')]);
 ui.block25([ui.h1('25% block width')]);
 ```

 ## FlexBox Grid
 Flexbox grid allow to divide a layout into multiple columns and rows. The Flexbox container take the full available width, and their height is determined by their inner content. A Flexbox layout has a direction in which child elements are laid out. The main axis is defined by rows or columns.

![Flexbox preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/flexbox.jpg)

 ```javascript
 ui.divH([ui.span('item 1'),ui.span('item 2'),ui.span('item 3')]); //rows
 ui.divV([ui.span('item 1'),ui.span('item 2'),ui.span('item 3')]); //collumns
 ui.divV([ui.span('item1'),ui.divH([ui.span('item 2'),ui.span('item 3')])]); //combined method with rows and columns
 ```

 ## Splitters
 Splitters - help to build the layout that contains several content areas. Each splitter contains the [box container](#boxes) which shrinks the content to a certain size.

 The splitters can specify by the horizontal or vertical orientation. In order to split vertically and horizontally at the same time, splitters need to be nested.

![Splitters preview](https://raw.githubusercontent.com/datagrok-ai/public/master/help/develop/img/splitters.jpg)

 ```javascript
 ui.splitH([ui.h1('Left'), ui.h1('Center'),  ui.h1('Right')])
 ui.splitV([ui.h1('Top'), ui.h1('Middle'), ui.h1('Bottom')])
 ui.splitH([ui.h1('Left'), ui.splitV([
   ui.h1('Right top'),
   ui.h1('Right bottom')
   ])
 ])
 ```

# Views
 ## Table View
 Table view - a view container with a grid table that contains a set of data that is structured in rows and columns. It allows the user to scroll in both directions and can contain large numbers of items and columns.
 ```javascript
 let table = grok.data.demo.demog();
 let view = grok.shell.addTableView(table);
 ```
 ## Simple View
 Simple view is an empty view container that can contain any kind of elements.
 ```javascript
 let view = grok.shell.newView('Simple View');
 ```
 ## Viewers
 ## Ribbon
 Ribbon is a layout container that appears in the header of the view.

 The ribbon panel can include text links and buttons, icons, dropdown menus, and Combobox, or any combination of those elements. Each collection with elements will be divided by a separator and evenly distributed in the container.
 ```javascript
 //Ribbon panel
 let view = DG.View.create();
 view.setRibbonPanels([
  [
  ui.divText('Custom panel')
  ],
  [
    ui.iconFA('search', () => grok.shell.info("clicked")),
    ui.iconFA('plus', () => grok.shell.info("plus"))
  ]
 ]);

 //Ribbon dropdown menu
 view.ribbonMenu = DG.Menu.create()
  .group('Menu')
  .item('element 1');
 ```

 ## Toolbox
 The toolbox is a container that appears in the left side of the view. It mostly used as a placement for UI controls and can contain any kind of elements such as buttons, fields, icons, links, dropdown menus, accordions, etc.
 ```javascript
 let view = grok.shell.newView('toolbox demo');
 let acc = ui.accordion();
 acc.addPane('header 1', () => ui.divText('Dynamic content'));
 view.toolbox = acc.root;
 ```

# Dialogs
The dialog is a control that informs about the task or contain some necessary information and requires decisions.
The dialog contains the following sections and options:
- Title: Title text appears in the dialog header.
- Content: This area contains the actual content of the dialog.
- Footer with actions: The footer can contain optional buttons and a context menu. If no buttons are defined, the default Close button is shown.

There are three types of dialogs:
- standard dialog
- modal dialog
- fullscreen modal dialog

For each dialog, you can set the position by viewport by x and y-asix.
```javascript
//standart dialog
ui.dialog('Standart dialog')
  .add(ui.span(['Some content...']))
  .onOK(() => {})
  .addContextAction('My Action', () => {}))
  .addButton('Optional button')
  .show({x: 300, y: 300});

//Fullscreen modal dialog
ui.dialog('Modal dialog')
  .add(ui.span(['Some content...']))
  .onOK(() => {})
  .showModal(true);

//modal dialog
ui.dialog('Modal dialog')
  .add(ui.span(['Some content...']))
  .onOK(() => {})
  .show();
```

# Elements
  ## Colors
  Color palette is a predefined set of colors. The colors are fixed and do not change.
  ```javascript
  let v = grok.shell.newView('palettes');

  function getBlock(c) {
    let block = ui.divText(DG.Color.toRgb(c));
    block.style.backgroundColor = DG.Color.toRgb(c);
    block.style.color = DG.Color.toRgb(DG.Color.getContrastColor(c));
    return block;
  }

  v.appendAll([
    ui.h1('Categorical palette with contrast text color'),
    ui.div(DG.Color.categoricalPalette.map(getBlock)),
    ui.h1('Category colors (looping over the palette)'),
    ui.div(DG.utils.identity(30).map(DG.Color.getCategoricalColor).map(getBlock))
  ]);
  ```
  ## Typography
  Typography sets default styles for headings, paragraphs, spans, and divs elements.
  ### Headers
  ```javascript
  ui.h1('Header 1');
  ui.h2('Header 2');
  ui.h3('Header 3');
  ```
  ### Paragraphs
  ```javascript
  ui.p('Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industrys standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.');
  ```
  ### Spans
  ```javascript
  ui.span(['span element']);
  ```
  ### Label
  A label is the name or title of a control or group of related controls.
  ```javascript
  ui.label('label text');
  ```
  ### Text blocks
  ```javascript
  ui.divText('Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industrys standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.');
  ```
  ## Tables
  Table - display information in a grid-like format of rows and columns. Table can contain interactive components such as buttons.
  ```javascript
  //Creates a visual table based on [items] and [renderer]
  var myList = [
  {key: 'first object', value: ui.button('button')},
  {key: 'second object', value: false},
  ];
  let table = DG.HtmlTable.create(myList, (item, idx) => [item.key, item.value]);

  //Creates a visual table based on [map]
  ui.tableFromMap({
    project: grok.shell.project.toMarkup(),
    time: new Date(),
  })
  ```
  ## Lists
  ```javascript
  ui.list([
  'element 1',
  grok.shell.user
  'element 3',
])
  ```
  ## Buttons
  Buttons allow users to trigger an action. There are 2 button types:
  - Regular button
  - Icon button
  - Big button (Use big button for make focus on main or single action)
  ```javascript
  ui.button('Regular button');
  ui.button(ui.iconFA('info'));
  ui.bigButton('Big button');
  ```
  ## Forms
  A form is used to present UI controls and allow to enter data in a structured way.
  ```javascript
  ui.inputs([
    ui.stringInput('Name'),
    ui.intInput('Age'),
    ui.buttonsInput([
      ui.bigButton('Apply'),
      ui.button('Cancel')
    ])  
  ]);
  ```
  ### Inputs
  Input field allows users to enter and edit text or numeric values in one line.
  There are two types of input fields:
  - **String input** allows to enter or edit text value
  - **Int input** allows to enter or edit numeric value
  ```javascript
  ui.stringInput('Name', 'Arthur Dent');
  ui.intInput('Age', 30);
  ```
  ### Text area
  The text area is an input control that allows the user to enter several lines of text.
  ```javascript
  ui.textArea('Text area text data');
  ```
  ### Dropdown selection
  The select control is used to select an item from a predefined list.
  ```javascript
  ui.choiceInput('Label', 'Value 1', ['Value 1', 'Value 2']);
  ```
  ### Selection
  The select control let's option to set a binary value (true/false). When the user clicks the selection control, it toggles between checked and unchecked.
  ```javascript
  ui.boolInput('Name', false);
  ```
  ### Group selection
  Group selection is commonly used to select one or more option from the predefined list.
  ```javascript
  ui.multiChoiceInput('Group label', ['Value 1', 'Value 2'], ['Value 3', 'Value 4']);
  ```
  ### Switch
  The toggle switch allows to set individual features to either active or inactive.
  ### Range slider
  Range slider is a UI control that enables to select a value range within a predefined interval.
  ```javascript
  ui.rangeSlider(0, 10, 2, 5);
  ```
  ## Icons
  The icon control displays the icon from the FontAwesome library. Icon can used as a button control. Use FontAwesome without 'fa' prefix.
  ```javascript
  ui.iconFA('question',()=>{})
  ```
  ## Images
# Components
  ## Accordions
  The accordion control is a container for grouping elements into panes. Each pane can collapse or expand and can contain any UI elements.
  ```javascript
  let acc = ui.accordion();
  acc.addPane('Pane label', () => ui.div([
    ui.h1('Header 1'),
    ui.button(' Button'),
    ui.textArea('multi\nline\ntext')])
  );
  ```

  ## Await (Loading indicator)
  Await informs the user about an ongoing operation.
  ```javascript
  ui.wait(async () => {
    let root = ui.div();
    return root;
  })
  ```

  ## Cards
  A card is used to display content composed of different elements whose size or supported actions can vary.
  ```javascript
  ui.card(
    ui.div([
    	ui.h1('header'),
      ui.divText('Some text'),
      ui.button('Regular')
    ]))
  ```

  ## Combo popup
  The combo box control allows users to select an item from a predefined list.
  ```javascript
  ui.div([ui.comboPopupItems('Combo popup label', {
        'Item 1': () => {},
        'Item 2': () => {},
    })]);
  ```

  ## Popup Menu (context menu)
  Menus appear upon interaction with a button, action, or other control. They display a list of choices, with one choice per lin. Menu can also have a multilevel list of choices.
  ```javascript
  let showMenu = () => {
    let showBalloon = (item) => grok.shell.info(item);
    DG.Menu.popup()
        .item("Show info", () => {}))
        .separator()
        .items(["First", "Second"], showBalloon)
        .show();
  };
  let text = ui.divText("Clickable element");
  text.addEventListener("click", showMenu);
  ```

  ## Property panel (Meta)
  Property panel is the right sidebar panel that used for showing the active item properties.
  ```javascript
  grok.shell.o = ui.h1('Property panel');
  ```

  ## Sidebar
  Sidebar control is a left side container that can contain items with toolbox pane.
  ```javascript
  grok.shell.sidebar.addPane('FIRST', () => ui.divText('A panel'), ui.iconFA('smile'));
  ```

  ## Tabs
  Tab control organize and allow navigation between a different content area or view.
  ```javascript
  ui.tabControl({
        'Fisrt Tab': ui.panel('First panel'),
        'Second Tab': () => ui.panel('Second panel')
    })
  ```

  ## Tag editor
  Tag editor is control that are small tag items that mainly serve to visualize selected items. Tags can be added, removed.
  ```javascript
  let editor = DG.TagEditor.create();
  editor.addTag('demo');
  editor.addTag('test');
  editor.addTag(1234);
  ```

  ## Taskbar progress
  Taskbar informs the user about loading progress at the bottom bar.
  ```javascript
  let pi = DG.TaskBarProgressIndicator.create('Progress...');
  setTimeout(() => {
    pi.close();
  }, 3000);
  ```

  ## Toasts
  A message toast is a small popup for success or error messages that disappears automatically after a few seconds.
  There are two types of message toast:
  - Success message
  - Error message
  ```javascript
  grok.shell.info('Success message')
  grok.shell.error('Error message');
  ```

  ## Top menu (Ribbon menu)
  ```javascript
  let view = grok.shell.newView('Demo View');
  view.ribbonMenu = DG.Menu.create()
    .group('Menu label')
    .item('Item label', () => {});
  ```  

  ## Tooltips
  Tooltip control display informative text when user hover over an element.
  ```javascript
  ui.tooltip.bind(ui.label('Label'), 'Tooltip message');
  ```

  ## Tree view
  A tree view control  presents a hierarchical view of information. Each item can have a number of subitems. One of the main use case is to display the hierarchically structured and to selecting one or more items out of a set of hierarchically structured items.
  ```javascript
  let tree = ui.tree();

  let group1 = tree.group('group 1', 1);
  group1.enableCheckBox(true);
  group1.item('item 1.1');
  group1.item('item 1.2');

  let subGroup1 = group1.group('group 1.1', 1.1, false);
  subGroup1.item('item 1.1.1');
  subGroup1.item('item 1.1.2');

  let group2 = tree.group('group 2', 2);
  group2.enableCheckBox();
  group2.item('item 2.1');
  group2.item('item 2.2').enableCheckBox();

  let groups = [group1, subGroup1, group2];
  let items = group1.items.slice()
    .concat(subGroup1.items)
    .concat(group2.items);
  ```
