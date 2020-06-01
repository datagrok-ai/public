// Creating a DataFrame from a JSON string

let table = DG.DataFrame.fromJson(`[
  {
    "name": "Roger Federer",
    "height": 185,
    "born": "August 8, 1981"
  },
  {
    "name": "Rafael Nadal",
    "height": 185,
    "born": "June 3, 1986"
  }
]`);

grok.shell.addTableView(table);
