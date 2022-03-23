import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';

import * as sql from 'sql.js';

//@ts-ignore: there are no types
import initSqlJs from '../sql-wasm.js';
import {_package} from '../package-test';
import {importSQLiteImpl} from '../import';

category('SQLite', () => {
  let chinookBytes: Uint8Array;
  let SQL: sql.SqlJsStatic;

  before(async () => {
    SQL = await initSqlJs({locateFile: () => _package.webRoot + 'dist/sql-wasm.wasm'});
    chinookBytes = await _package.files.readAsBytes('chinook.sqlite');
  });

  test('import', async () => {
    const chinookTables = importSQLiteImpl(chinookBytes, SQL).map((t) => t!.name);
    const expectedTables = [
      'albums', 'sqlite_sequence', 'artists', 'customers', 'employees', 'genres', 'invoices',
      'invoice_items', 'media_types', 'playlists', 'playlist_track', 'tracks', 'sqlite_stat1',
    ];

    expect(chinookTables.length, expectedTables.length);
    expect(chinookTables.filter((tableName) => !expectedTables.includes(tableName)).length, 0);
  });
});
