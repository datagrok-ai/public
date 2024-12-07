let t = grok.shell.table('Med Chem') ?? grok.data.demo.molecules();
let rnd = (n) => Math.floor(n * Math.random());
let rnd1 = (n) => Math.ceil(n * Math.random());
let rndCat = (categories) => categories[rnd(categories.length)];
const addCatColumn = (columnName, categories, tags) => {
  let c = t.columns.addNewString(columnName).init((i) => rndCat(categories));
  if (tags != null) {
    for (let key of Object.getOwnPropertyNames(tags))
      c.tags[key] = tags[key];
  }
  return c;
};

// HTS image
addCatColumn('HTS', ['https://upload.wikimedia.org/wikipedia/commons/thumb/4/43/Housebeemouth100x.jpg/120px-Housebeemouth100x.jpg', 'https://en.wikipedia.org/wiki/File:RiceStemcs400x1.jpg', 'https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/Rabbitttestis100x2.jpg/120px-Rabbitttestis100x2.jpg', 'https://upload.wikimedia.org/wikipedia/commons/thumb/5/52/FernProthallium400x.jpg/120px-FernProthallium400x.jpg']);

// HTS data
t.columns.addNew('HTS data', DG.TYPE.DATA_FRAME).init((i) => grok.data.demo.doseResponse(96));

let stages = t.columns.addNewString('Stage').init((i) => 'Stage ' + rnd1(4));
let species = addCatColumn('Species', ['Mouse', 'Rat']);
let strains = {
  'Mouse': ['C57BL/6', 'BALB/c', 'Swiss Webster', 'DBA/2'],
  'Rat': ['Sprague Dawley', 'Wistar', 'Long-Evans', 'Fisher 344', 'Lewis'],
};
t.columns.addNewString('Strain').init((i) => rndCat(strains[species.get(i)]));

// Chemistry
const chemTags = {quality: 'Molecule', 'cell.renderer': 'Molecule'};
addCatColumn('Core', ['n1nc([*:3])n([*:2])c1[*:1]', 'C1C2CN([*:2])CC2CN1[*:1]', 'CC1(N[*:1])CCN([*:2])C1'], chemTags);
addCatColumn('R1', ['N#Cc1cccc(C[*:1])c1', 'COCC(=O)[*:1]', 'CC(C)N1CCC(N2CCN([*:1])CC2)CC1', 'O=C(C1CC1)[*:1]', 'C#CCN1CCN([*:1])CC1', 'COCCOC(C)[*:1]', 'CCNC(=O)CN1CCN([*:1])CC1', 'COCCOC1CCN([*:1])C1', 'Cc1cnc(C(=O)[*:1])cn1', 'CC(=O)[*:1]'], chemTags);
addCatColumn('R2', ['N#CC[*:2]', 'O=C(c1ccon1)[*:2]', 'Cn1nnnc1C[*:2]', 'O=C(C[*:2])NCC(F)(F)F', 'CS(=O)(=O)CC(=O)[*:2]', 'O=C(C[*:2])NCC(F)(F)F', 'O=C(C[*:2])NCC(F)(F)F', 'c1csc(C[*:2])n1', 'c1csc(C[*:2])n1', 'c1nonc1C[*:2]', 'O=C(CCOCC(F)(F)F)[*:2]'], chemTags);
addCatColumn('R3', ['C1CCC([*:3])OC1', 'Brc1ccc([*:3])s1', 'O=C(Nc1ccccn1)[*:3]', 'c1ccc(C[*:3])nc1', 'c1ccc(-n2cc([*:3])cn2)cc1', 'c1ccc([*:3])cc1', 'CN(Cc1ccc(S(C)(=O)=O)cc1)[*:3]', 'c1cncc([*:3])c1', 'CON(C)[*:3]', 'Cc1ccc(SC[*:3])cc1', 'Brc1ccc([*:3])s1', 'Cc1cccc([*:3])n1', 'C1CC1C[*:3]', 'C1CC([*:3])C1', 'CC1(C)CCN([*:3])CC1', 'CN(Cc1nc2ccccc2s1)[*:3]', 'CC(O)(C1CCN([*:3])CC1)C(F)(F)F'], chemTags);

// Classifications
addCatColumn('BSEP classification', ['>= 60 uM', '<= 30 uM', 'Inconclusive']);
addCatColumn('HLM CLint classification', ['<=25 uL/min/mg', '>100 uL/min/mg', 'Inconclusive']);
addCatColumn('RLM CLint classification', ['<=30 uL/min/mg', '>120 uL/min/mg', 'Inconclusive']);
addCatColumn('LE-MDCK Classification', ['<=5 cm-6/s', '>=5 cm-6/s', 'Inconclusive']);
addCatColumn('PAMPA Classification', ['<=-5.3 cm/s', '>= -4.5 cm/s', 'Inconclusive']);
addCatColumn('pH6.8 HT Solubility Classification', ['<10 uM', '>100 uM', 'Inconclusive']);


// Random boolean columns
const boolNames = ['hERG alert', 'Patent pending', 'Completed'];
for (let name of boolNames)
  t.columns.addNewBool(name).init((i) => rnd(3) === 1);


const warnings = ['Mutagenicity', 'Tumorigenicity', 'Irritative effects', 'Reproductive effects'];
const warningsCol = t.columns.addNewString('Toxicity').init((i) => warnings.filter((i, j) => rnd(3) === 1).join(','));
warningsCol.setTag(DG.TAGS.CELL_RENDERER, 'Tags');

t.columns.addNewString('Idea name').init((i) => 'GRK-' + rnd(10000));
addCatColumn('Idea status', ['Created', 'Approved', 'Synthesized']);
addCatColumn('Idea author', ['Jack Hammer', 'Paul Smith', 'Maria DePiotra', 'Marlon Voicek']);
let ideaCreated = t.columns.addNewDateTime('Idea created');
let ideaApproved = t.columns.addNewDateTime('Idea approved');
let ideaSynthesized = t.columns.addNewDateTime('Idea synthesized');
for (let i = 0; i < t.rowCount; i++) {
  ideaCreated.set(i, dayjs().year(1990 + rnd(30)).month(rnd1(12)).day(rnd1(28)).startOf('day'));
  ideaApproved.set(i, ideaCreated.get(i).add(rnd(100), 'day'));
  ideaSynthesized.set(i, ideaApproved.get(i).add(rnd(100), 'day'));
}
ideaCreated.tags.format = 'yyyy/MM/dd';
ideaApproved.tags.format = 'yyyy/MM/dd';
ideaSynthesized.tags.format = 'yyyy/MM/dd';

t.meta.detectSemanticTypes();
grok.shell.addTableView(t);