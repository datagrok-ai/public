export const layoutConf = {
  innerRadius: 250,
  outerRadius: 300,
  cornerRadius: 0,
  gap: 0.04,
  labels: {
    display: true,
    position: 'center',
    size: '14px',
    color: '#000000',
    radialOffset: 30,
  },
  ticks: {
    display: false,
    color: 'grey',
    spacing: 10000000,
    labels: true,
    labelSpacing: 10,
    labelSuffix: '',
    labelDenominator: 1000000,
    labelDisplay0: false,
    labelSize: '10px',
    labelColor: '#000000',
    labelFont: 'default',
    majorSpacing: 5,
    size: {
      minor: 2,
      major: 5,
    }
  },
  events: {}
};

export function topSort(graph) {
  let stack = [];

  function visit(id, node, stack) {
    node.visited = true;
    for (let target of node.targets) {
      if (!graph[target].visited) visit(target, graph[target], stack);
    }
    stack.unshift(node.datum);
  }

  for (let [id, node] of Object.entries(graph)) {
    if (!node.visited) visit(id, node, stack);
  }

  return stack;
}
