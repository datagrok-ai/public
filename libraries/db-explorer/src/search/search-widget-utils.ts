import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/** processes a table view that is going to be inside of power search result */
export function processPowerSearchTableView(tv: DG.TableView) {
  const maxDepth = 6;
  let curRoot: HTMLElement | null = tv.root;
  for (let i = 0; i < maxDepth; i++) {
    if (curRoot == null)
      break;
    if (curRoot.classList.contains('power-pack-widget-host'))
      break;

    curRoot = curRoot.parentElement;
  }

  if (curRoot == null || !curRoot.classList.contains('power-pack-widget-host'))
    return;

  curRoot.style.setProperty('height', 'calc(100vh - 200px)', 'important');
  const header: HTMLElement | null = curRoot.querySelector('.d4-dialog-header');
  if (!header)
    return;

  const iEl: HTMLElement | null = header.querySelector('i');
  if (!iEl)
    return;
  const addTVButton = ui.icons.add(() => {
    const layout = tv.saveLayout();
    const addedTV = grok.shell.addTableView(tv.dataFrame);
    addedTV.loadLayout(layout);
    setTimeout(() => grok.shell.v = addedTV);
  }, 'Add To Workspace');
  addTVButton.style.width = '18px';
  header.insertBefore(addTVButton, iEl);
}

export type PowerSearchTQOptions = {
  layout?: DG.ViewLayout
  postProcess?: (tv: DG.TableView) => void | Promise<void>
};

export function powerSearchQueryTable(queryCall: DG.FuncCall, options?: PowerSearchTQOptions): DG.Widget {
  return DG.Widget.fromRoot(ui.wait(async () => {
    await queryCall.call();
    const df: DG.DataFrame = queryCall.getOutputParamValue();
    if (!df || (df.rowCount ?? 0) === 0)
      return ui.divText('No results found');

    await df.meta.detectSemanticTypes();
    await grok.data.detectSemanticTypes(df);
    const tv = DG.TableView.create(df, false);
    setTimeout(async () => {
      processPowerSearchTableView(tv);
      tv._onAdded();
      if (options?.postProcess)
        await options.postProcess(tv);
      if (options?.layout)
        tv.loadLayout(options.layout);
    }, 200);
    return tv.root;
  }));
}


// takes a pattern which contains some constants and variable regions (marked with ${})
// tries to parse it and extract the values of the variables from the query
export function matchAndParseQuery(matchPattern: string, query: string): string[] | undefined | null {
  //matchPattern is a string like "First 1000 rows of ${} from ${}"
  try {
    query ??= '';
    query = query.trim();
    const result: string[] = [];
    // replace the ${} in matchPattern with the regex that matches any string
    const regexPattern = matchPattern.replace(/\$\{\}/g, '(.*)');
    // make it case insensitive
    const regex = new RegExp(regexPattern, 'i');
    const match = query.match(regex);
    if (!match)
      return null;
    for (let i = 1; i < match.length; i++) {
      if (!match[i])
        return null;
      result.push(match[i]);
    }
    return result;
  } catch (_e) {
    return undefined;
  }
}
