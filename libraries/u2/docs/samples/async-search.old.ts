import * as ui from 'datagrok-api/ui';

interface Compound {
  id: string;
  name: string;
  mw: number;
}

export function compoundSearch(search: (query: string) => Promise<Compound[]>): HTMLElement {
  const results = ui.div();
  let latest = 0;
  let timer: ReturnType<typeof setTimeout> | undefined;
  const load = async (query: string) => {
    const call = ++latest;
    ui.setUpdateIndicator(results, true);
    try {
      const items = await search(query);
      if (call !== latest)
        return;
      ui.empty(results);
      results.append(items.length === 0 ? ui.divText('No compounds match') :
        ui.table(items, (c) => [c.id, c.name, c.mw.toFixed(1)], ['Id', 'Name', 'MW']));
    }
    catch (e) {
      if (call !== latest)
        return;
      ui.empty(results);
      results.append(ui.divText(e instanceof Error ? e.message : String(e)), ui.button('Retry', () => load(query)));
    }
    finally {
      if (call === latest)
        ui.setUpdateIndicator(results, false);
    }
  };
  const query = ui.input.search('Search', {onValueChanged: (q) => {
    clearTimeout(timer);
    timer = setTimeout(() => load(q), 150);
  }});
  load('');
  return ui.divV([ui.form([query]), results]);
}
