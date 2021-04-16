// https://github.com/jasondavies/newick.js
export function parseNewick(a) {for(var e=[],r={},s=a.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){var n=s[t];switch(n){case"(":var c={};r.branchset=[c],e.push(r),r=c;break;case",":var child={};e[e.length-1].branchset.push(child),r=child;break;case")":r=e.pop();break;case":":break;default:var h=s[t-1];")"==h||"("==h||","==h?r.name=n:":"==h&&(r.length=parseFloat(n))}}return r};

export function newickToDf(newick, filename) {
  let parent = 'root';
  let i = 0;
  const obj = parseNewick(newick);
  const nodes = [], parents = [], distances = [];

  function traverse(obj) {
    if (obj === null || typeof obj != "object" ) return;
    if (!Array.isArray(obj)) {
      let name = obj.name;
      if (!name) name = `node-${i}`, i += 1;
      nodes.push(name);
      distances.push(obj.length);
      parents.push(parent);
      parent = name;
    }
    Object.values(obj).forEach(value => traverse(value));
  }
  traverse(obj);

  const df = DG.DataFrame.fromColumns([
    DG.Column.fromList('string', 'node', nodes),
    DG.Column.fromList('string', 'parent', parents),
    DG.Column.fromList('double', 'distance', distances),
  ]);
  df.name = `df-${filename.slice(0, -4)}`;
  return df;
};
