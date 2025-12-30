const fruits = [{name: 'apple', price: 5}, {name: 'banana', price: 6}, {name: 'orange', price: 7}];

grok.shell.topMenu
  .item('Foo!', () => grok.shell.info('Foo clicked'))
  .separator()
  .group('Custom')
  .items(fruits, item => grok.shell.info(item), {
    getTooltip: item => `${item.name} ${item.price}`,
    isChecked: item => item.name === 'apple',
    isValid: item => item.name === 'orange' ? 'Oranges are not allowed' : null,
    onMouseEnter: item => grok.shell.info(item.name),
    toString: item => item.name
  })
  .endGroup();
