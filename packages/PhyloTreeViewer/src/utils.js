export function newickToDf(newick, filename) {
  let parent = null;
  let i = 0;
  const parsedNewick = d3.layout.newick_parser(newick);
  const obj = parsedNewick.json;
  const nodes = [], parents = [], distances = [], annotations = [];

  function traverse(obj) {
    if (obj === null || typeof obj != 'object' ) return;

    let name = obj.name;
    if (!name) name = obj.name = `node-${i}`, i++;

    nodes.push(name);
    distances.push(obj.attribute ? parseFloat(obj.attribute) : null);
    annotations.push(obj.annotation);
    parents.push(parent);

    if (!obj.children) return;
    const childrenNum = obj.children.length;
    const prevParent = parent;
    parent = name;

    for (let i = 0; i < childrenNum; i++) {
      traverse(obj.children[i]);
      if (i === childrenNum - 1) parent = prevParent;
    }
  }
  traverse(obj);

  const columns = [
    DG.Column.fromList('string', 'node', nodes),
    DG.Column.fromList('string', 'parent', parents),
  ];

  if (distances.some(d => d !== null)) {
    columns.push(DG.Column.fromList('double', 'distance', distances));
  }

  if (annotations.some(a => !!a)) {
    columns.push(DG.Column.fromList('string', 'annotation', annotations));
  }

  const df = DG.DataFrame.fromColumns(columns);

  df.name = `df-${filename.slice(0, -4)}`;
  df.setTag('.newick', newick);
  df.setTag('.newickJson', JSON.stringify(parsedNewick));

  return df;
};

// https://stackoverflow.com/questions/5525071/how-to-wait-until-an-element-exists
export function waitForElm(id, checkFrequency = 100, timeout = 1000) {
  const startTime = Date.now();
  return new Promise((resolve, reject) => {
    (function loopSearch() {
      const element = document.getElementById(id);

      if (element) {
        return resolve(element);
      } else {
        setTimeout(() => {
          if ((Date.now() - startTime) > timeout) {
            reject('Timeout');
          }
          loopSearch();
        }, checkFrequency);
      }
    })();
  });
}
