import * as DG from 'datagrok-api/dg';

import {compareIssues, GitHubIssue, isOpen} from './issue';

const GREEN = 0xFFE6F4EA;
const AMBER = 0xFFFFF3E0;
const RED = 0xFFFDE7E7;
const ORANGE = 0xFFFCE7D6;
const YELLOW = 0xFFFCF7DA;
const GREY = 0xFFF3F3F3;

const COLUMN_COLORS: Record<string, Record<string, number>> = {
  status: {'Done': GREEN, 'Blocked': RED, 'In Progress': AMBER, 'In Review': AMBER, 'Testing': AMBER,
    'Ready for QA': AMBER},
  priority: {Blocker: RED, Urgent: RED, High: ORANGE, Medium: YELLOW},
  resolution: {'Done': GREEN, 'Unresolved': AMBER},
};

function issuesToDataFrame(issues: GitHubIssue[], name: string): DG.DataFrame {
  const sorted = [...issues].sort(compareIssues);
  const col = (colName: string, get: (i: GitHubIssue) => string) => DG.Column.fromStrings(colName, sorted.map(get));
  const df = DG.DataFrame.fromColumns([
    col('ticket', (i) => `<a href="${i.url}" target="_blank">${i.key}</a>`),
    col('summary', (i) => i.title),
    col('status', (i) => i.boardStatus ?? (isOpen(i) ? 'Open' : 'Closed')),
    col('priority', (i) => i.priority ?? ''),
    col('assignee', (i) => i.assignees.join(', ')),
    col('resolution', (i) => i.resolution),
    col('type', (i) => i.type ?? ''),
    col('client', (i) => i.clients.join(', ')),
    col('milestone', (i) => i.milestone ?? ''),
  ]);
  df.name = name;
  df.getCol('ticket').setTag(DG.TAGS.CELL_RENDERER, 'html');
  return df;
}

export function issuesGrid(issues: GitHubIssue[], name: string): DG.Grid {
  const df = issuesToDataFrame(issues, name);
  const grid = DG.Viewer.grid(df, {showColumnGridlines: false, allowBlockSelection: false});
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell)
      return;
    const colors = COLUMN_COLORS[gc.gridColumn.name];
    if (colors)
      gc.style.backColor = colors[gc.cell.valueString] ?? GREY;
  });
  grid.root.classList.add('ua-issues-grid');
  const ticket = grid.col('ticket');
  if (ticket)
    ticket.width = 90;
  const summary = grid.col('summary');
  if (summary)
    summary.width = 360;
  return grid;
}
