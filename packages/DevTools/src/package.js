import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

async function loadSnippets() {
  return (await grok.dapi.scripts.list({ filter: '#demo and (language in ("javascript"))' }));
}

function getSnippets(ent, snippets) {
  return snippets.slice(0, 3);
}

//tags: autostart
export async function describeCurrentObj() {
  const apiSnippets = await loadSnippets();

  grok.events.onAccordionConstructed.subscribe(acc => {
    const ent = acc.context;

    if (ent) {
      const snippets = getSnippets(ent, apiSnippets);
      const snippetNames = snippets.map(s => ui.divText(s.name, { classes: 'd4-link-action' }));
      let editor = ui.textInput('', '');
      editor.input.style = 'width: 0; height: 0; visibility: hidden;';
      editor.root.style.display = 'none';

      snippetNames.forEach((el, idx) => el.addEventListener('click', () => {
        editor.value = snippets[idx].script;
        editor.root.style.display = 'block';
        editor.input.style = 'width: 200; height: 300; visibility: visible;';
      }));

      let snippetsPane = acc.getPane('Snippets');      
      if (!snippetsPane) snippetsPane = acc.addPane('Snippets', () => {
        return ui.divV([
          ui.divText(ent.name),
          ...snippetNames,
          editor.root
        ]);
      });
    }
  });
}
