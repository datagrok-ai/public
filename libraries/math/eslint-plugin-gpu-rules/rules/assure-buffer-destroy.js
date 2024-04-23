module.exports = {
  meta: {
    type: 'problem',
    docs: {
      description: 'Ensure all created GPU buffers are destroyed',
      category: 'Best Practices',
      recommended: false
    },
    schema: [] // No options for this rule
  },
  create(context) {
    const buffers = new Map();

    return {
      CallExpression(node) {
        // console.log(node.parent.id.name); //bufferName
        // console.log(node.callee.object.name) //bufferName in case of destroy;
        // console.log(node.callee.object.name); //device
        // console.log(node.callee.property && node.callee.property.name); //createBuffer/destroy()
        // Check for buffer creation
        if (node.parent && node.parent.id && node.parent.id.name && node.callee &&
          node.callee.property && node.callee.property.name === 'createBuffer' &&
          node.callee.object && node.callee.object.name === 'device'
        )
          buffers.set(node.parent.id.name, node);

        // Check for buffer destruction
        if (node.callee && node.callee.property && node.callee.property.name === 'destroy' &&
          node.callee.object && node.callee.object.name
        )
          buffers.delete(node.callee.object.name);
      },
      onCodePathEnd: function(codePath, node) {
        if (codePath.origin && codePath.origin == 'program') {
          if (buffers.size > 0) {
            buffers.forEach((reportNode, bufferName) => {
              context.report({
                node: reportNode,
                message: `Buffer ${bufferName} is not destroyed`
              });
            });
          }
        }
      },

    };
  }
};
