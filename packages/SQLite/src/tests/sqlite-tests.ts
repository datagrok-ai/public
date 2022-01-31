import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {importSQLite} from '../package';
import {symmetricDifference} from 'set-operations';

export const _package = new DG.Package();

category('SQLite', () => {
  let chinookBytes: Uint8Array;

  before(async () => {
    chinookBytes = await _package.files.readAsBytes('chinook.sqlite');
  });

  test('import', async () => {
    const chinookTables = (await importSQLite(chinookBytes)).map((t) => t!.name);
    const expectedTables = [
      'albums', 'sqlite_sequence', 'artists', 'customers', 'employees', 'genres', 'invoices',
      'invoice_items', 'media_types', 'playlists', 'playlist_track', 'tracks', 'sqlite_stat1',
    ];

    expect(symmetricDifference(chinookTables, expectedTables, true).length, 0);
  });
});
