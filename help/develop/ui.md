<!-- TITLE: Grok UI -->
<!-- SUBTITLE: -->

<!-- This is a user-centric view on the Datagrok applications UI development-->

# Layouts

 ## Starting point
 To start building a layout you need either a [View](#simple-view) (in most cases) or a [Dialog](#dialogs).
 Every View or Dialog is a [panel container](#panels)

 ## Containers
This is a simple container. It can contain any elements, such as inputs, images, text, etc.
It doesn't have own height. The container height depends on its children.
```javascript
ui.div([ui.h1('Header'), ui.p('Paragraph text'), 'just text', DG.Viewer.scatterPlot(grok.data.demo.demog())])
```
If you place a container within a [box container](#boxes) , it will inherit the box size and show scroll bars automatically.
 ## Boxes
This is a fixed-size container. It doesn't depend on children element sizes, but shrinks them to certain size.

 ## Panels
 The panel is a simple container similar to the [Containers](#Containers). It has full available wide and its height depends on its children. Also, panels have 10px margins for all sides.
 ```javascript
 ui.panel([ui.h1('Header'), ui.p('Paragraph text'), 'just text', DG.Viewer.scatterPlot(grok.data.demo.demog())])
 ```

 ## Blocks

 ## Splitters
 Splitters - help to build the layout that contains several content areas. Each splitter contains the [box container](#boxes) which shrinks the content to a certain size.

 The splitters can specify by the horizontal or vertical orientation. In order to split vertically and horizontally at the same time, splitters need to be nested.

 ```javascript
 ui.splitH([ui.h1('Left'), ui.h1('Center'),  ui.h1('Right')])
 ui.splitV([ui.h1('Top'), ui.h1('Middle'), ui.h1('Bottom')])
 ui.splitH([ui.h1('Left'), ui.splitV([ui.h1('Right top'),ui.h1('Right bottom')])])
 ```

# Views
 ## Table View
 ## Simple View
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

# Dialogs

# Elements
  ## Colors
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
  ### Text blocks
  ```javascript
  ui.divText('Lorem Ipsum is simply dummy text of the printing and typesetting industry. Lorem Ipsum has been the industrys standard dummy text ever since the 1500s, when an unknown printer took a galley of type and scrambled it to make a type specimen book.');
  ```
  ## Tables
  ## Lists
  ## Buttons
  ## Forms
  ### Form input
  ### Form selection
  ### Form group selection
  ### Form switch
  ### Form range slider
  ## Icons
  ## Images
# Components
  ## Accordions
  ## Cards
  ## Menu
  ## Top menu
  ## Popup menu
  ## Tabs
  ## Tag editor
  ## Sidebar
  ## Toasts
  ## Tooltips
