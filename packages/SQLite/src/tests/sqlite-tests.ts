import {before, category, expect, test} from '@datagrok-libraries/test/src/test';

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

  // The bundled SQLite 3.36 rejected the whole schema of a file with a RIGHT JOIN view, so no table opened
  test('import file with RIGHT JOIN view', async () => {
    const source = new SQL.Database();
    source.run(`CREATE TABLE a(id INTEGER, x TEXT); CREATE TABLE b(id INTEGER, y TEXT);
      INSERT INTO a VALUES (1, 'x1'); INSERT INTO b VALUES (1, 'y1'), (2, 'y2');
      CREATE VIEW v_data_daily AS SELECT a.id, a.x, b.y FROM a RIGHT JOIN b ON a.id = b.id;`);
    const bytes = source.export();
    source.close();

    const tables = importSQLiteImpl(bytes, SQL);
    expect(tables.map((t) => t.name).join(','), 'a,b');
    expect(tables[1].rowCount, 2);
  });
});
