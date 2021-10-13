/**
 * RDKit-based substructure filters that uses Datagrok's collaborative filtering.
 *
 * 1. On onRowsFiltering event, only FILTER OUT rows that do not satisfy this filter's criteria
 * 2. Call dataFrame.rows.requestFilter when filtering criteria changes.
 * */

class SubstructureFilter extends DG.Filter {

  constructor() {
    
    super();
    this.WHITE_MOL = `
  0  0  0  0  0  0  0  0  0  0999 V2000
M  END
`;
    this.molfile = '';
    this.root = ui.divV([]);
    this.root.appendChild(grok.chem.sketcher((_, molfile) => {
        this.molfile = molfile;
        this.dataFrame.rows.requestFilter();
      }));
      
  }

  attach(dFrame) {
    
    this.dataFrame = dFrame;
    this.column = this.dataFrame.columns.bySemType(DG.SEMTYPE.MOLECULE);
    this.dataFrame.onRowsFiltering.subscribe((_) => this.applyFilter());
  
  }

  reset() {
    
    this.dataFrame.filter.setAll(true, false);
    if (this.column?.tags['chem-scaffold-filter']) {
      delete this.column.tags['chem-scaffold-filter'];
    }
    this.dataFrame.filter.fireChanged();
  
  }

  applyFilter() {
    
    if (!this.molfile || this.molfile.endsWith(this.WHITE_MOL)) {
      this.reset();
      return;
    }
    grok.functions.call('Chem:searchSubstructure', {
      'molStringsColumn' : this.column,
      'molString': this.molfile,
      'substructLibrary': false,
      'molStringSmarts': ''
    })
    .then((bitset_col) => {
      this.dataFrame.filter.copyFrom(bitset_col.get(0));
      this.column.setTag('chem-scaffold-filter', this.molfile);
      this.dataFrame.filter.fireChanged();
    }).catch((e) => {
      console.warn(e);
      this.reset();
    })
    
  }
}