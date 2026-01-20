/**
 * Atom manipulations
 * Original example: http://www.cheminfo.org/Chemistry/Cheminformatics/OpenChemLib_js/index.html
 * For more details, visit https://github.com/cheminfo/openchemlib-js
 */

let v = grok.shell.newView('atoms');
let smiles = 'c1(ccc2N=C(C)N(C(=O)c2c1)c3ccc(OC)cc3)NC(=S)Nc4ccccc4';
let m = OCL.Molecule.fromSmiles(smiles);
let result = {};
let atoms = [];

result.atoms = atoms;
for (let i = 0; i < m.getAllAtoms(); i++) {
  let atom = {};
  atoms.push(atom);
  atom.abnormalValence = m.getAtomAbnormalValence(i);
  atom.charge = m.getAtomCharge(i);
  atom.radical = m.getAtomRadical(i);
  atom.cipParity = m.getAtomCIPParity(i);
  atom.color = m.getAtomColor(i);
  atom.esrGroup = m.getAtomESRGroup(i);
  atom.esrType = m.getAtomESRType(i);
  atom.atomicNo = m.getAtomicNo(i);
  atom.customLabel = m.getAtomCustomLabel(i);
  atom.label = m.getAtomLabel(i);
  atom.list = m.getAtomList(i);
  atom.listString = m.getAtomListString(i);
  atom.mapNo = m.getAtomMapNo(i);
  atom.mass = m.getAtomMass(i);
  atom.parity = m.getAtomParity(i);
  atom.queryFeatures = m.getAtomQueryFeatures(i);
  atom.radical = m.getAtomRadical(i);
  atom.x = m.getAtomX(i);
  atom.y = m.getAtomY(i);
  atom.z = m.getAtomZ(i);
  atom.allConnected = m.getAllConnAtoms(i);
  atom.hydrogens = m.getAllHydrogens(i);
  atom.implicitHydrogens = m.getImplicitHydrogens(i);
  atom.connected = m.getConnAtoms(i);
  atom.pi = m.getAtomPi(i);
  atom.ringSize = m.getAtomRingSize(i);
  atom.occupiedValence = m.getOccupiedValence(i);
  atom.freeValence = m.getFreeValence(i);
  atom.stabilized = m.isStabilizedAtom(i);
  atom.ringBondCount = m.getAtomRingBondCount(i);
  atom.stereoBond = m.getStereoBond(i);

  atom.isConfigurationUnknown = m.isAtomConfigurationUnknown(i);
  atom.isParityPseudo = m.isAtomParityPseudo(i);
  atom.isStereoCenter = m.isAtomStereoCenter(i);
  atom.isNaturalAbundance = m.isNaturalAbundance(i);
  atom.isPurelyOrganic = m.isPurelyOrganic(i);
  atom.isAutoMapped = m.isAutoMappedAtom(i);
  atom.isAllylic = m.isAllylicAtom(i);
  atom.isAromatic = m.isAromaticAtom(i);
  atom.isRing = m.isRingAtom(i);
  atom.isFlatNitrogen = m.isFlatNitrogen(i);
  atom.isSimpleHydrogen = m.isSimpleHydrogen(i);
  atom.isSmallRing = m.isSmallRingAtom(i);
}

let bonds = [];
result.bonds = bonds;
for (let i = 0; i < m.getAllBonds(); i++) {
  let bond = {};
  bonds.push(bond);
  bond.ringSize = m.getBondRingSize(i);
  bond.cipParity = m.getBondCIPParity(i);
  bond.esrGroup = m.getBondESRGroup(i);
  bond.esrType = m.getBondESRType(i);
  bond.length = m.getBondLength(i);
  bond.order = m.getBondOrder(i);
  bond.parity = m.getBondParity(i);
  bond.queryFeatures = m.getBondQueryFeatures(i);
  bond.bridgeMinSize = m.getBondBridgeMinSize(i);
  bond.bridgeMaxSize = m.getBondBridgeMaxSize(i);
  bond.type = m.getBondType(i);
  bond.typeSimple = m.getBondTypeSimple(i);
  bond.isBridge = m.isBondBridge(i);
  bond.isParityPseudo = m.isBondParityPseudo(i);
  bond.isAromaticBond = m.isAromaticBond(i);
  bond.isParityUnknownOrNone = m.isBondParityUnknownOrNone(i);
  bond.isStereo = m.isStereoBond(i);
  bond.atom1 = m.getBondAtom(0, i);
  bond.atom2 = m.getBondAtom(1, i);
}

v.append(grok.chem.svgMol(smiles));

grok.shell.addTableView(DG.DataFrame.fromObjects(atoms));
grok.shell.addTableView(DG.DataFrame.fromObjects(bonds));