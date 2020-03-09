// HTML-based grid cells with the markup-driven content.
// For code-driven rendering, check out "html-dynamic-cell.js" sample.
// See also: https://datagrok.ai/help/overview/markup

let view = grok.addTableView(grok.testData('demog', 5000));
let col = view.grid.columns.byName('disease');
view.grid.options({'rowHeight': 100});
col.width = 200;
col.format =
    `<div style="display:flex; flex-direction:column; padding: 5px">
  <div>
    <div>sex: <b>#{t.row[sex]}</b></div>
    <div> height: <span style="background-color:#{t.color(height)}"><b>#{t.row[height]}</b></span> </div>
    <button class="btn d4-btn-ok btn-raised">CONTACT</button>
  </div>
</div>`;
col.cellType = 'html';