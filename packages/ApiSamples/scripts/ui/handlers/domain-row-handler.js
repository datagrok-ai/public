// Custom ObjectHandler for domain-table rows (DG.DomainRow).
//
// A domain row's semantic type is its entity type name, '<schema>.<table>'.
// Registering an ObjectHandler with that type overrides the platform's generic
// domain-row rendering for that table only (see the Grit package's grit.issue
// handler). DG.DomainRow exposes: schemaName, tableName, typeName, semValue,
// values (raw column map), id, version, authorId, createdOn, updatedOn.

const ISSUE_TYPE = 'grit.issue';

class DomainIssueHandler extends DG.ObjectHandler {
  get type() {return ISSUE_TYPE;}

  // Wins dispatch for exactly the rows of this table.
  isApplicable(x) {
    return (x instanceof DG.DomainRow && x.typeName === ISSUE_TYPE) ||
      (x instanceof DG.SemanticValue && x.semType === ISSUE_TYPE);
  }

  renderCard(x) {
    const row = x instanceof DG.SemanticValue ? x.value : x;
    const status = ui.span([row.values.status ?? '']);
    status.style.color = '#1f8fff';
    return ui.divV([
      ui.divH([ui.divText(row.semValue), status]),
      ui.divText(row.values.title ?? ''),
    ], 'd4-gallery-item');
  }
}

// Register: wins forEntity dispatch (context panel, cards, entity view).
DG.ObjectHandler.register(new DomainIssueHandler());

// Discover handlers registered for a semType (e.g. from the Grit package).
DG.ObjectHandler.forSemType(ISSUE_TYPE).then((handlers) =>
  grok.shell.info(`${handlers.length} handler(s) for ${ISSUE_TYPE}`));
