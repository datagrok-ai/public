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
 The panel is a simple container similar to the [Containers]. It has full available wide and its height depends on its children. Also, panels have 10px margins for all sides.
 ```javascript
 ui.panel([ui.h1('Header'), ui.p('Paragraph text'), 'just text', DG.Viewer.scatterPlot(grok.data.demo.demog())])
 ```

 ## Blocks
 ## Splitters

# Views
 ## Table View
 ## Simple View
 ## Viewers
 ## Ribbon
 ## Toolbox

# Dialogs

# Content
 ## Typography
  ### Headers
  ### Paragraphs
  ### Spans
 ## Inputs
 ## Forms
 ## Icons
 ## Buttons
 ## Links
 ## Colors
 ## Images
 ## Tab Control
 ## Menu
 ## Cards
 ## Tables
