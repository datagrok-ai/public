function toHierarchy(dataFrame, splitByColumnNames) {
  let data = {
    name: 'All',
    children: []
  };

  let aggregated = dataFrame.groupBy(splitByColumnNames).aggregate();
  let columns = aggregated.columns.byNames(splitByColumnNames);
  let parentNodes = columns.map(c => null);

  for (let i = 0; i < aggregated.rowCount; i++) {
    let idx = i === 0 ? 0 : columns.findIndex((col) => col.get(i) !== col.get(i - 1));

    for (let colIdx = idx; colIdx < columns.length; colIdx++) {
      let parentNode = colIdx === 0 ? data : parentNodes[colIdx - 1];
      let node = {name: columns[colIdx].getString(i)};
      parentNodes[colIdx] = node;
      if (!parentNode.children)
        parentNode.children = [];
      parentNode.children.push(node);
    }
  }

  return data;
};

grok.shell.info(toHierarchy(grok.data.demo.demog(), ['sex', 'race']));