//api: DG.U2.component
// A function component: component() runs the builder with the new control's scope as `self`; answering
// {root, ...api} gives the control that API (getters stay live), and `name` is its automation id.
const counter = DG.U2.component((self) => {
  const count = DG.U2.signal(0);
  const text = ui.divText('0');
  self.effect(() => text.textContent = `clicked ${count.value} times`);
  return {root: ui.divV([text, ui.button('Click', () => count.value++)]), get count() {
    return count.value;
  }};
}, {name: 'counter'});
grok.shell.newView('u2 component', [counter.root]);
grok.shell.info(`${counter.root.dataset.u2Name}: ${counter.count}`);
