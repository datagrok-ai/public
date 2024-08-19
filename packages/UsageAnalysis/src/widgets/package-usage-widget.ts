import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {getUsersTable} from '../tabs/packages';

export class PackageUsageWidget extends DG.Widget {
  date: string = 'this week';
  groupsP: Promise<DG.Group[]> = grok.dapi.groups.getGroupsLookup('All users');

  constructor(pack: DG.Package) {
    super(ui.box(null, {classes: 'ua-widget ua-package-widget'}));
    this.root.append(ui.waitBox(async () => {
      const df: DG.DataFrame = await this.groupsP.then((groups) => grok.functions.call('UsageAnalysis:PackagesUsage',
        {date: this.date, groups: [groups[0].id], packages: ['all']}));
      df.rows.removeWhere((r: DG.Row) => r.get('package') !== pack.name);
      const usersHistogram = df
        .groupBy(['uid'])
        .sum('count')
        .aggregate();
      if (usersHistogram.rowCount === 0)
        return ui.divText('No data to display', {style: {color: 'var(--failure)'}});
      return getUsersTable(usersHistogram);
    }));
    this.root.append(ui.label(this.date));
  }
}
