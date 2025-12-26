//wait, waitBox, loader

async function sleep(time) {
  return new Promise(resolve => setTimeout(resolve, time));
}

//load indicator
let loader = ui.loader();

//wait container
let wait = ui.wait(async () => {
  await sleep(500);
  return ui.divText('ui.wait element');
});

//wait container. Recommend to use inside tabControl, accordion, box, view.box
let waitBox = ui.waitBox(async () => {
  await sleep(1000);
  return ui.divText('ui.waitBox element');
});

grok.shell.newView('Demo, loader, wait, waitBox', [wait, waitBox, loader]);