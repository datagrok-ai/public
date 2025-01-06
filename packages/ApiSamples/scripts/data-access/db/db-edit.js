// Real-time updates in the DB as you edit the table

const dbEdit = (t) => {
  const dbTableName = t.tags['DbTable'];
  const dbIdColumnName = 'id';

  let v = grok.shell.addTableView(t);
  v.grid.onCellValueEdited.subscribe(async (c) => {
    const idColumnName = 'id'
    const dbColumnName = c.cell.column.tags['DbColumn'];
    const sql =
      'update ' + dbTableName +
      ' set ' + dbColumnName + " = '" + c.cell.value + "'" +
      ' where ' + dbIdColumnName + " = '" + c.cell.row.get(idColumnName) + "'";

    grok.shell.info('Executing ' + sql);
 
    await grok.data.db.query('System:Datagrok', sql);
    grok.shell.info('Success!'); 
  });
}

let result = await grok.data.db.query('System:Datagrok', 'select * from models');
dbEdit(result);