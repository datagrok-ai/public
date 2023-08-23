//wait, waitBox, loader

async function sleep() {
  return new Promise(resolve => setTimeout(resolve, 2000));
}

//load indicator
let loader = ui.loader();

//wait container
let wait = ui.wait(async () => {
  await sleep();
  return ui.divText('ui.wait element');
})

//wait container. Recommend to use inside tabControl, accordion, box, view.box
let waitBox = ui.waitBox(async () => {
  await sleep();
  return ui.divText('ui.waitBox element');
})

wait.style.border = '1px solid var(--grey-1)';
waitBox.style.border = '1px solid var(--grey-1)';
grok.shell.newView('Demo, loader, wait, waitBox', [wait, waitBox, loader]);