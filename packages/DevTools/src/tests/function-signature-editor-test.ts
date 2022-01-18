import * as DG from "datagrok-api/dg";
import * as grok from "datagrok-api/grok";
import * as ui from "datagrok-api/ui";
import {category, before, test, after, expect, delay} from "@datagrok-libraries/utils/src/test";

const delayForRender = 500;

category('FSE exists', ()=> {
  before(async ()=> {
    const v = DG.View.createByType('ScriptView');
    grok.shell.addView(v);
    await delay(delayForRender);
  })

  test('FSE button exists', async () => {
    const fseButton = grok.shell.v.getRibbonPanels().flat().find(
      (elem) => !!elem.querySelector('.fa-magic')
    ).lastChild as HTMLElement;
    if (!fseButton) throw "Failed: FSE button does not exist"
  })

  after(async () => {
    grok.shell.v.close();
  })
})

category('FSE opens', () => {
  before(async ()=> {
    const v = DG.View.createByType('ScriptView');
    grok.shell.addView(v);
    await delay(delayForRender);
  })

  test('FSE button flips view', async () => {
    const fseButton = grok.shell.v.getRibbonPanels().flat().find(
      (elem) => !!elem.querySelector('.fa-magic')
    ).lastChild as HTMLElement;
    fseButton.click();

    await delay(delayForRender);
    const currentView = grok.shell.v;
    const allTabsExist = 
      !!currentView.root.querySelector('[name="PROPERTIES"]') &&
      !!currentView.root.querySelector('[name="PARAMETERS"]') &&
      !!currentView.root.querySelector('[name="CODE"]') &&
      !!currentView.root.querySelector('[name="UI"]') 
    expect(allTabsExist, true);
  })

  after(async () => {
    grok.shell.v.close();
  })
});

category('FSE closes', () => {
  before(async ()=> {
    const v = DG.View.createByType('ScriptView');
    grok.shell.addView(v);
    await delay(delayForRender);
  })

  test('FSE button flips view once more', async () => {
    const fseButton = grok.shell.v.getRibbonPanels().flat().find(
      (elem) => !!elem.querySelector('.fa-magic')
    ).lastChild as HTMLElement;
    fseButton.click();
    await delay(delayForRender);

    const currentView = grok.shell.v;
    const editorButton = currentView.getRibbonPanels().flat().find(
      (elem) => !!elem.querySelector('.fa-code')
    ).lastChild as HTMLElement;
    editorButton.click();
    await delay(delayForRender);

    const editorExists = !!currentView.root.querySelector('.CodeMirror.cm-s-default');
    expect(editorExists, true);
  })

  after(async () => {
    grok.shell.v.close();
  })
})