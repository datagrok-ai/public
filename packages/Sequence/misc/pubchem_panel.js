// class _ProteinViewerInfoPanelProvider extends InfoPanelProvider<Cell> {
//   String get name => 'PDB Molecule';
//   String get description => 'Visualizes a molecule by its PDB id';
//
//   static RegExp uniprotProtein = new RegExp('http://purl.uniprot.org/uniprot/(.*)');
//   static int _id = 0;
//
//   Future<bool> isApplicable(dynamic x) async => x is Cell && x.isCell && x.column.type == Types.STRING && x.value != null &&
// (x.column.name == 'pdb' || uniprotProtein.hasMatch(x.value));
//
//   static Future<String> uniprotToPdb(String uniprotUri) async {
//   var sparql = '''PREFIX up:<http://purl.uniprot.org/core/>
//   PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
//   SELECT ?protein ?db
//   WHERE
// {
// ?protein a up:Protein .
// ?protein rdfs:seeAlso ?db .
// ?db up:database <http://purl.uniprot.org/database/PDB>.
//   filter (?protein = <$uniprotUri>)
//   }''';
//
//   var result = await Sparql.query(Sparql.UNIPROT, sparql);
//
//   if (result.rowCount == 1 && result['db'] != null)
//   return result['db'][0];
//
//   return null;
//   }
//
//   Future<InfoPanel> makePanel(Cell cell, Map<String, dynamic> params) async {
//     var m = uniprotProtein.firstMatch(cell.value);
//     var pdb = m == null ? cell.value : await uniprotToPdb(cell.value);
//     var viewer = new NglViewer('ngl_viewer$_id', init: (NglViewer v) => v.load(url: 'rcsb://$pdb'));
//
//     return new InfoPanel('PDB Molecule', viewer.root);
//   }
//     }