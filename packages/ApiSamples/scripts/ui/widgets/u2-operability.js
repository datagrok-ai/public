//api: DG.U2.operability
// What a user is offered but cannot use: DG.U2.operability lists each interactive element that is covered, takes
// no pointer, sits in an inert region with no reason, is disabled with no reason, is named one thing and tooltipped
// another, or is an options-rail child taller than the rail.
const dead = ui.button('Edit', () => grok.shell.info('edited'));
dead.style.pointerEvents = 'none';
const off = ui.button('Save', () => {});
off.setAttribute('aria-disabled', 'true');
const why = ui.button('Submit', () => {});
why.setAttribute('aria-disabled', 'true');
why.setAttribute('aria-description', 'Name the campaign first');
const root = ui.divV([dead, off, why]);
grok.shell.newView('u2 operability', [root]);

for (const f of DG.U2.operability(root, {hitTest: true}))
  grok.shell.warning(`${f.rule}: ${f.detail}`);
