// https://datagrok.ai/help/visualize/viewers/markup

let user = grok.shell.user;

let table = grok.data.testData('demog', 10000);
let view = grok.shell.addTableView(table);

let sp = DG.Viewer.scatterPlot(table, {
  xColumnName: 'age',
  yColumnName: 'weight',
  markerDefaultSize: 5,
  markerType: 'circle',
  colorColumnName: 'site',
  legendVisibility: 'Never'
});

let hist = DG.Viewer.histogram(table, {
  valueColumnName: 'weight'
});

let boxPlot = DG.Viewer.boxPlot(table, {
  valueColumnName: 'weight',
  categoryColumnNames: ['race']
});

let filter = DG.Viewer.filters(table, {
  columnNames: ['race']
});

let markup = `<div style="padding:0 20px">
<h2>HTML Markup</h2>
Markup Viewer lets you combine HTML markup with the properties that are dynamically evaluated against
the current dataset. This is useful for creating custom dashboards, or telling stories with data. Here are some
examples that illustrate most important concepts:

    <h3>HTML elements</h3>
Use the standard HTML markup to add elements, or <i>style the document</i>.
<img src="/images/ribbon/project.svg"/>
    <br>
    Custom objects: <span>${user.toMarkup()}</span>

<h3>Metadata</h3>
<p> Access metadata, including names and tags for tables and columns.</p>
Name: <span>#{t.name}</span> <br>
Tags: Description: <span>#{t.tags[${DG.TAGS.DESCRIPTION}]}</span> <br>

<h3> Counts </h3>
<p> Try filtering out rows to see that the values are synchronized: </p>
Row count: <span>#{t.rowCount}</span> <br>
Selected: <span>#{t.selection.trueCount}</span> <br>
Filter: <span>#{t.rows.filters}</span> <br>
Filtered: <span>#{t.filter.trueCount}</span> <br>
<br>

<h3>Referencing data:</h3>
Current row: <span>#{t.currentRow}</span> <br>
Current column: <span>#{t.currentCol}</span> <br>
Current cell: <span>#{t.currentCell}</span> <br>
Current subject's race: <span>#{t.currentRow[race]}</span> <br>
Formula evaluation: current subject's BMI = weight / (height^2) = <span>#{formula(\${weight} / (\${height} *\${height} / 10000 ))}</span> <br>
Editable age: <span>#{t.editor[age]}</span> <br>

<h3>Descriptive statistics:</h3>
avg(age): <span>#{t.stats.avg(age)}</span> <br>
stdev(age): <span>#{t.stats.stdev(age)}</span> <br> <br>

<h3>Interactivity</h3>
<button class="btn" data-onclick="post http://mycompany/report/#{t.currentRow[subj]}">
    Report protocol violation for <span>#{t.currentRow[subj]}</span> </button>
calls external REST API, parametherized by the subject id

</div>`;

view.markup({content: markup});
