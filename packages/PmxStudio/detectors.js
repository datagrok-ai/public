/**
 * The class contains semantic type detectors.
 * Detectors are functions tagged with `DG.FUNC_TYPES.SEM_TYPE_DETECTOR`.
 * See also: https://datagrok.ai/help/develop/how-to/define-semantic-type-detectors
 * The class name is comprised of <PackageName> and the `PackageDetectors` suffix.
 * Follow this naming convention to ensure that your detectors are properly loaded.
 */
class PmxstudioPackageDetectors extends DG.Package {
  // A study table (its columns' roles in the 'pmx.study' tag, set by Monitor and Analyze) gets PMX Studio's items in
  // its grid's context menu, and its table views PMX Studio's painting and menu. This file loads with the platform:
  // it only reads the tag and builds the items; the package loads when an item is chosen or a study's view opens.

  //name: autostart
  //meta.role: autostart
  //description: PMX Studio's items in a study table's context menu, and its marks in the study's table views
  autostart() {
    const roles = (df) => {
      try {
        const r = JSON.parse((df && df.getTag('pmx.study')) || 'null');
        return r && df.col(r.animal) && df.col(r.day) ? r : null;
      } catch (e) {
        return null;
      }
    };

    grok.events.onContextMenu.subscribe((e) => {
      const cell = e.args.item;
      if (!(cell instanceof DG.GridCell) || !cell.isTableCell || cell.tableRowIndex == null)
        return;
      const df = cell.grid.dataFrame;
      const r = roles(df);
      if (!r)
        return;
      const i = cell.tableRowIndex;
      const animal = String(df.get(r.animal, i));
      const selected = df.selection.get(i) && df.selection.trueCount > 1 ?
        Array.from(df.selection.getSelectedIndexes()) : null;
      const these = selected || [i];
      const one = selected ? `${selected.length} selected measurements` : `${animal} on day ${df.get(r.day, i)}`;
      const ofAnimal = [];
      for (let k = 0; k < df.rowCount; k++) {
        if (String(df.get(r.animal, k)) === animal)
          ofAnimal.push(k);
      }
      // the marks, columns of the table (study-table.ts): a tick, or a file's own text (true, yes, 1, x)
      const cols = [df.col('Excluded'), df.col('Flagged')].filter((c) => c);
      const on = (v) => v === true || v === 1 || (typeof v === 'string' && /^(true|yes|y|1|x)$/i.test(v.trim()));
      const marked = (rows) => rows.some((k) => cols.some((c) => on(c.get(k))));
      const act = (rows, action, what) => () => grok.functions.call('PmxStudio:pmxMarkStudyRows',
        {table: df, rows: rows.join(','), action, what});
      const m = e.args.menu.group('PMX Studio');
      m.item(`Flag ${one}…`, act(these, 'flag', one));
      m.item(`Exclude ${one}…`, act(these, 'exclude', one));
      m.item(`Flag animal ${animal}…`, act(ofAnimal, 'flag', `animal ${animal}`));
      m.item(`Exclude animal ${animal}…`, act(ofAnimal, 'exclude', `animal ${animal}`));
      if (marked(these))
        m.item(`Clear the mark on ${one}`, act(these, 'clear', one));
      if (marked(ofAnimal))
        m.item(`Clear the marks on animal ${animal}`, act(ofAnimal, 'clear', `animal ${animal}`));
    });

    grok.events.onViewAdded.subscribe((view) => {
      if (view instanceof DG.TableView && roles(view.dataFrame))
        grok.functions.call('PmxStudio:pmxStudyTableView', {view});
    });
  }
}
