/** A message's plural forms, by `Intl.PluralRules` category: `other` always, the others where the language has them. */
export type PluralForms = {
    readonly other: string;
} & {
    readonly [K in Exclude<Intl.LDMLPluralRule, 'other'>]?: string;
};
/** A message: text, a template with `{name}` placeholders, or plural forms whose `{count}` is the number. */
export type Message = string | PluralForms;
/**
 * What a host gives a sketcher to translate it (I18N-002): any of the catalog's messages by name, in the catalog's form,
 * and element names as `element.<Symbol>` (`element.N`: `"Stickstoff"`). Left out, a message is English.
 */
export type SketcherStrings = Readonly<Record<string, Message>>;
/** The Mod key as a tooltip names it (key names never translate): Ctrl, or ⌘ on a Mac. */
export declare const KEY_NAMES: Readonly<{
    ctrl: "Ctrl";
    cmd: "⌘";
}>;
export declare const STRINGS: {
    /** The canvas's accessible name. */
    readonly canvasLabel: 'Molecule drawing';
    /** The toolbars' groups in the default configuration (A11Y-003). */
    readonly groupHistory: 'History';
    readonly groupCanvas: 'Canvas';
    readonly groupBonds: 'Bonds';
    readonly groupElements: 'Elements';
    /** The ring bar's group (spike rings-chains, reading 7; Ketcher's has no label). */
    readonly groupRings: 'Rings';
    /** The left bar's group of the hand tool and the eraser, as Ketcher's first card holds them (spike view, V4). */
    readonly groupHandErase: 'Hand and eraser';
    /** The left bar's first group since spike select (S1): the selection palette, and its tools, Ketcher's names. */
    readonly groupSelection: 'Selection';
    readonly selectionTools: 'Selection tools';
    readonly rectangleSelection: 'Rectangle selection';
    readonly lassoSelection: 'Lasso selection';
    readonly structureSelection: 'Structure selection';
    /** The top bar's zoom group (spike view, UI-001). */
    readonly groupZoom: 'Zoom';
    /** The top bar's group of commands on the whole structure, Ketcher's chemistry commands (spike layout, UI-001). */
    readonly groupStructure: 'Structure';
    /** The top bar's clipboard group (spike clipboard, C3). */
    readonly groupClipboard: 'Clipboard';
    readonly stereoBonds: 'Stereo bonds';
    /** The controls (A11Y-002, reading 14). */
    readonly undo: 'Undo';
    readonly redo: 'Redo';
    readonly clearCanvas: 'Clear canvas';
    readonly eraser: 'Eraser';
    readonly singleBond: 'Single bond';
    readonly doubleBond: 'Double bond';
    readonly tripleBond: 'Triple bond';
    readonly aromaticBond: 'Aromatic bond';
    readonly wedgeBond: 'Wedge bond';
    readonly hashedWedgeBond: 'Hashed wedge bond';
    readonly wavyBond: 'Wavy bond';
    readonly crossedDoubleBond: 'Crossed double bond';
    readonly periodicTable: 'Periodic table';
    /** The ring tools, named after their rings (spike rings-chains, A11Y-002), and the chain tool (reading 23). */
    readonly benzene: 'Benzene';
    readonly cyclopentadiene: 'Cyclopentadiene';
    readonly cyclohexane: 'Cyclohexane';
    readonly cyclopentane: 'Cyclopentane';
    readonly cyclopropane: 'Cyclopropane';
    readonly cyclobutane: 'Cyclobutane';
    readonly cycloheptane: 'Cycloheptane';
    readonly cyclooctane: 'Cyclooctane';
    readonly chain: 'Chain';
    /** The R-group tools (spike rgroups, RG7): the left bar's group and palette, and their tools, Ketcher's names. */
    readonly groupRGroups: 'R-groups';
    readonly rgroupTools: 'R-group tools';
    readonly rgroupLabel: 'R-group label';
    readonly attachmentPoint: 'Attachment point';
    /** The S-Group tool (spike sgroups, G3): the left bar's group after the R-groups, and its tool, Ketcher's name. */
    readonly groupSGroups: 'S-groups';
    readonly sgroupTool: 'S-Group';
    /** The charge tools (spike charges-checks, K4; CHG-001): the left bar's group after the R-groups, and its two tools. */
    readonly groupCharges: 'Charges';
    readonly chargePlus: 'Charge plus';
    readonly chargeMinus: 'Charge minus';
    /** The Structure Library (spike templates, RING-013): the ring bar's last button and its dialog (Ketcher's title). */
    readonly structureLibrary: 'Structure Library';
    /** The Structure Library while its module loads, and one that did not load (milestone-4-fixes, F10). */
    readonly libraryLoading: 'Loading the Structure Library…';
    readonly libraryNotLoaded: 'The Structure Library could not be loaded';
    /** Its search field (RING-015, reading 9), and what it shows when nothing matches (reading 10). */
    readonly search: 'Search';
    readonly searchHint: 'Name, group or abbreviation';
    readonly noItemsFound: 'No items found';
    /** A group's header: its name and how many templates it shows ("Aromatics (18)"). */
    readonly groupCount: '{group} ({count})';
    /**
     * Spike abbreviation-library (L1 to L3): the library's tabs, as Ketcher names them, and the tab list's name. Ketcher's own
     * names of its lists, which the lookup's options say each option is (L6).
     */
    readonly libraryTabs: 'Libraries';
    readonly templateLibrary: 'Template Library';
    readonly functionalGroups: 'Functional Groups';
    readonly saltsAndSolvents: 'Salts and Solvents';
    /** The lookup (L6, ABBR-005, ABBR-007; Ketcher's words) and Choose Nickname (L5, ABBR-004; ChemDraw's). */
    readonly abbreviationLookup: 'Abbreviation lookup';
    readonly chooseNickname: 'Choose Nickname';
    readonly lookupHint: 'Type an element, a template or a group';
    readonly noMatchingResults: 'No matching results';
    readonly lookupElement: 'Element';
    readonly lookupTemplate: 'Template';
    readonly lookupGroup: 'Functional group';
    readonly lookupSalt: 'Salt or solvent';
    /** The hand tool (V4, Ketcher's) and the zoom group's controls (spike view, A11Y-002). */
    readonly hand: 'Hand';
    readonly zoomIn: 'Zoom in';
    readonly zoomOut: 'Zoom out';
    readonly zoomReset: 'Zoom to 100%';
    readonly zoomFit: 'Fit to view';
    readonly zoomLevel: 'Zoom level';
    /**
     * The zoom as one control (spike bar-fit, B7, Ketcher's): its menu's name, the button's name with the level it shows, and
     * the menu's name for the zoom to 100% (Ketcher's "Actual size").
     */
    readonly zoomMenu: 'Zoom';
    readonly zoomButton: 'Zoom: {percent}';
    readonly zoomActualSize: 'Actual size';
    /** Layout and Clean Up (spike layout, A11Y-002): Ketcher's names. */
    readonly layout: 'Layout';
    readonly cleanUp: 'Clean Up';
    /**
     * What the status line says when the engine cannot lay out a selection around the atoms it holds (spike layout-tools,
     * L2): the depictor refused (two pieces of one ring system held), or no such layout keeps the stereo.
     */
    readonly layoutRefused: 'This selection cannot be laid out around the rest: select whole rings, or the whole structure';
    readonly layoutLosesStereo: 'No layout of this selection around the rest keeps its stereochemistry';
    /**
     * Layout's and Clean Up's tooltips (spike layout-tools, L9: each acts on the selection when there is one, else on
     * everything), which say how they differ (spike bar-fit, B8, the product owner's: layout.md L2, L3).
     */
    readonly layoutTip: 'Layout: redraw from scratch and fit the view (the selection, if any)';
    readonly cleanUpTip: 'Clean Up: tidy bonds and angles, keeping the drawing in place (the selection, if any)';
    /** Clean Up's menu (B8): its name, its opener's too, and the line under each item. */
    readonly cleanUpMenu: 'Clean Up and Layout';
    readonly cleanUpLine: 'Tidy bonds and angles, keep the drawing in place';
    readonly layoutLine: 'Redraw from scratch, fit the view';
    /** Straighten (spike layout-tools, L5, L9): RDKit's straightenDepiction, in the selection's and the canvas's menus (P1: no button, no key). */
    readonly straighten: 'Straighten';
    /** Cut, Copy and Paste (spike clipboard, C3, A11Y-002). */
    readonly cut: 'Cut';
    readonly copy: 'Copy';
    readonly paste: 'Paste';
    /** The zoom as the top bar shows it. */
    readonly percent: '{n}%';
    /** Between two keys of one tooltip: "Ctrl+Delete or Ctrl+Backspace". */
    readonly keysOr: ' or ';
    /** A tooltip's hotkey, in brackets after the name (UI-007). The toolbars set the key in a `<kbd>` (fix round 4). */
    readonly withHotkey: '{name} ({key})';
    /** A disabled tool's tooltip: no tool is disabled since spike stereo, which enabled the wavy and crossed bonds (BOND-001). */
    readonly notAvailableYet: '{name}: not available yet';
    /** The inline label editor (ATOM-012). */
    readonly labelEditor: 'Atom label';
    /**
     * Why an edit is refused (basics-fixes F3), in the canvas's polite live region. An atom on the canvas is named by its
     * element, as the status line speaks of the atom under the pointer ("this C", spike charges-checks, K14), and one the
     * gesture adds by its element alone (a template's attachment atom, "Cl"); one message for each thing there is no room for.
     */
    readonly noRoomBond: 'No room on {atom} for another bond';
    readonly noRoomRing: 'No room on {atom} for the ring';
    readonly noRoomChain: 'No room on {atom} for the chain';
    readonly noRoomTemplate: 'No room on {atom} for the template';
    /**
     * A template the edit operations cannot make (spike reports-oct5, R4): why it is not placed, when picked and at every
     * gesture. Its dative bonds are no reason since spike engine-dative-ez (D2): the edit operations make them.
     */
    readonly templateRadicals: 'Not placed: this template has more radical electrons on an atom than the sketcher can draw';
    readonly templateBondType: 'Not placed: this template has a bond of a type the sketcher cannot draw';
    readonly noRoomMerge: 'No room on {atom} for the merge';
    /** An atom on the canvas, as a refusal names it: "this C" (K14). */
    readonly thisAtom: 'this {symbol}';
    /**
     * A label's radical mark on an atom with no hydrogen to give for it (spike charges-checks, K15, provisional): "C." on a
     * carbon with four bonds. Named by the element the label makes.
     */
    readonly noRoomForRadical: 'No room on this {symbol} for a radical';
    /** A bond the Hand would drop within half a bond of another atom (M3). */
    readonly tooClose: 'Too close to {atom}';
    /**
     * A problem of the molecule (spike charges-checks, K10; CHECK-002), in ChemDraw's words: the status line says it while
     * the pointer or the keyboard's hotspot is on a flagged atom, and `warnings` gives it as each warning's `message`. An
     * over-valent unlabelled carbon; any other atom with an invalid valence, by its element; an aromatic system with no
     * Kekulé structure; a label kept as an alias (ATOM-016, CHECK-005).
     */
    readonly problemCarbon: 'Too many bonds to this unlabeled carbon';
    readonly problemValence: 'Invalid valence for this {symbol}';
    readonly problemAromatic: 'These aromatic bonds have no Kekulé structure';
    readonly problemParentheses: "Parentheses don't match";
    readonly problemLabel: "Crux Sketch can't interpret this label";
    /** A stereocentre nothing specifies in a structure drawn with stereo bonds (spike stereo, S13; CHECK-007). */
    readonly problemStereo: 'Unspecified stereocentre';
    /** An abbreviation with more bonds than attachments (spike abbreviations, A2, K1; P-row, provisional), which RDKit draws whole. */
    readonly problemAbbreviation: 'Too many bonds to this abbreviation';
    /**
     * The formula and mass readout (spike charges-checks, K12; INFO-001 to INFO-009, UI-030): its name and its toggle's,
     * the caption over a selection's values, the masses' labels, what it says of a query (ChemDraw's), and the note beside
     * the masses when R labels or attachment points count in the formula but not in the masses (INFO-006).
     */
    readonly infoPanel: 'Formula and mass';
    readonly infoSelection: 'Selection';
    readonly infoWeight: 'MW';
    readonly infoWeightName: 'Molecular weight';
    readonly infoExactMass: 'Exact mass';
    readonly infoQuery: 'Formula cannot be computed for queries';
    readonly infoRNotCounted: 'R not counted';
    /**
     * A flip or a 180° rotation of a part attached to the rest of its structure by more than one bond (spike rotate-flip,
     * R4): refused, as F3 says a refusal.
     */
    readonly cannotFlip: {
        readonly one: 'Cannot flip a part attached by {count} bond';
        readonly other: 'Cannot flip a part attached by {count} bonds';
    };
    readonly cannotTurn: {
        readonly one: 'Cannot rotate by 180° a part attached by {count} bond';
        readonly other: 'Cannot rotate by 180° a part attached by {count} bonds';
    };
    /** `.` over an atom that has both attachment points already (spike keyboard-a11y, K4, KEY-015). */
    readonly noFreeAttachment: 'No free attachment point on {atom}';
    /**
     * What a paste says in the status line when it changes nothing (spike clipboard, C2, C3; CLIP-012): text that is no
     * structure; text in a format a paste does not read yet, by its name (InChI, CDXML, SMARTS); a clipboard the browser
     * does not let the buttons use, with the key that works instead (CLIP-015).
     */
    readonly pasteNothing: 'Nothing to paste: the clipboard holds no structure';
    readonly pasteNotYet: 'Pasting {format} is not supported yet';
    readonly clipboardBlocked: '{action} is not permitted here: press {key} instead';
    /**
     * Copy As and Paste Special (spike copy-as, P1, P2): the menus the Copy and Paste buttons open, by ChemDraw's names. A
     * refused write, and text that is not what Paste Special was told to read.
     */
    readonly copyAs: 'Copy As';
    readonly pasteSpecial: 'Paste Special';
    readonly copyAsBlocked: 'Copy As {format} is not permitted here';
    /**
     * What a Copy As did not write (milestone-4-fixes, F4), as ChemDraw warns on Copy As: "Copied as SMARTS; not written:
     * enhanced stereo groups"; a MOL V2000 written as V3000, and why: "Copied as MOL V3000: MOL V2000 cannot hold enhanced
     * stereo groups"; the losses by the loss table's codes.
     */
    readonly copiedLost: 'Copied as {format}; not written: {what}';
    readonly copiedV3000: 'Copied as MOL V3000: MOL V2000 cannot hold {why}';
    readonly copiedV3000Lost: 'Copied as MOL V3000: MOL V2000 cannot hold {why}; not written: {what}';
    readonly copiedV3000Needed: 'Copied as MOL V3000, which this drawing needs';
    readonly lossAlias: 'aliases';
    readonly lossRLabelNumbers: 'combined R labels';
    readonly lossAttachmentPoints: 'attachment points';
    readonly lossStereoGroups: 'enhanced stereo groups';
    readonly lossChiralFlag: 'the chiral flag';
    readonly lossQuery: 'query features';
    readonly lossSgroup: 'S-groups';
    readonly lossHydrogenBond: 'hydrogen bonds';
    readonly lossZeroOrder: 'zero-order bonds';
    readonly lossAtomMaps: 'atom maps';
    readonly lossStereo: 'some stereo';
    readonly lossRLabel: 'R labels';
    readonly lossDummy: 'dummy atoms and pseudoatoms';
    readonly lossDative: 'dative bonds';
    readonly lossOther: 'some features';
    /** Copy As of a part a format cannot hold at all (spike inchi, N6): nothing is copied, and the status line says why. */
    readonly copiedNone: 'Not copied: {format} cannot hold {what}';
    readonly pasteNotAs: 'Nothing to paste: the clipboard holds no {format}';
    readonly pasteNotAsMol: 'Nothing to paste: the clipboard holds no molfile';
    /** The Paste Special dialog (P2), where the browser does not let the sketcher read the clipboard. */
    readonly pasteSpecialHint: 'The browser does not let the sketcher read the clipboard. Paste the text here with {key}, choose how to read it, then OK.';
    readonly pasteText: 'Text to paste';
    readonly readAs: 'Read as';
    readonly notReadAs: 'This text cannot be read as {format}';
    readonly notReadAsMol: 'This text cannot be read as a molfile';
    /** A file dropped on the canvas that adds nothing (P4, CLIP-019): no structure in it, or a format not read yet. */
    readonly dropNothing: 'Nothing to add: {name} holds no structure';
    readonly unnamedFile: 'the file';
    readonly dropNotYet: 'Dropping {format} files is not supported yet';
    /** A reaction pasted or dropped (spike formats, F12, IO-031): refused, saying so; a written one is the engine's error, in its words. */
    readonly reactionsRefused: 'Reactions are not supported';
    /**
     * ChemDraw's text (spike cdxml, X3, X5): a dropped binary CDX file, which is not read; a paste with no text, which is what
     * a page sees of ChemDraw's own copy, with the way to copy CDXML text there; the CDXML module that did not load.
     */
    readonly cdxNotRead: 'CDX files are not read: save as CDXML';
    readonly pasteNoText: 'Nothing to paste as a structure; in ChemDraw use Copy As > CDXML Text ({key})';
    readonly cdxmlNotLoaded: 'CDXML cannot be read or written here: the CDXML module did not load';
    /** A CDXML paste or drop while the CDXML module loads (milestone-4-fixes, F2): said until the read report replaces it. */
    readonly cdxmlLoading: 'Loading the CDXML reader…';
    /** Copy As CDXML while the CDXML module loads (milestone-5-fixes; the ui-critic's milestone-5 note): said until it is written. */
    readonly cdxmlWriterLoading: 'Loading the CDXML writer…';
    /**
     * InChI (spike inchi, N2, N4): a paste or drop while the InChI module loads (as CDXML's, F2); the module that did not load;
     * an InChIKey pasted or dropped, which no structure can be made of.
     */
    readonly inchiLoading: 'Loading the InChI reader…';
    /** Copy As InChI or InChIKey while the InChI module loads (milestone-5-fixes; the ui-critic's milestone-5 note): said until it is written. */
    readonly inchiWriterLoading: 'Loading the InChI writer…';
    readonly inchiNotLoaded: 'InChI cannot be read or written here: the InChI module did not load';
    readonly inchiKeyRefused: 'An InChIKey cannot be drawn: it is a hash, not a structure';
    /**
     * What a CDXML read left out, in the status line (spike cdxml, X3, the read report): "Pasted 2 structures; left out: 1
     * text, 1 arrow", after a paste, a drop or a host's write; the labels kept as aliases after it.
     */
    readonly readPasted: {
        readonly one: 'Pasted {count} structure';
        readonly other: 'Pasted {count} structures';
    };
    readonly readDropped: {
        readonly one: 'Added {count} structure';
        readonly other: 'Added {count} structures';
    };
    readonly readWritten: {
        readonly one: 'Read {count} structure';
        readonly other: 'Read {count} structures';
    };
    readonly readLeftOut: '{read}; left out: {what}';
    readonly readAlso: '{report}; {more}';
    readonly leftStructure: {
        readonly one: '{count} structure that could not be read';
        readonly other: '{count} structures that could not be read';
    };
    readonly leftText: {
        readonly one: '{count} text';
        readonly other: '{count} texts';
    };
    readonly leftGraphic: {
        readonly one: '{count} graphic';
        readonly other: '{count} graphics';
    };
    readonly leftArrow: {
        readonly one: '{count} arrow';
        readonly other: '{count} arrows';
    };
    readonly leftPicture: {
        readonly one: '{count} picture';
        readonly other: '{count} pictures';
    };
    readonly leftBracket: {
        readonly one: '{count} bracket';
        readonly other: '{count} brackets';
    };
    readonly keptAliases: {
        readonly one: '{count} label kept as an alias';
        readonly other: '{count} labels kept as aliases';
    };
    /**
     * What any read dropped, named in the read report by the engine's code (milestone-4-fixes, F1): "Added 1 structure; left
     * out: the records after the first".
     */
    readonly lostRecords: 'the records after the first';
    readonly lostSdData: 'SD data fields';
    readonly lostRgroupDefinitions: 'R-group definitions';
    readonly lostCoordinates3d: '3D coordinates (laid out in 2D)';
    readonly lostSgroup: 'S-groups it cannot keep';
    readonly lostQuery: 'query features it cannot keep';
    readonly lostCxSection: 'CXSMILES sections it cannot keep';
    readonly lostMolfile: 'molfile properties it cannot keep';
    readonly lostOther: 'some of the text';
    /** The buttons floating beside a selection (spike rotate-flip, R5, XFORM-021): their group, and each one's name. */
    readonly selectionActions: 'Selection actions';
    readonly flipHorizontal: 'Flip horizontal';
    readonly flipVertical: 'Flip vertical';
    readonly deleteSelection: 'Delete selection';
    /** The angle a rotation has turned, beside its handle (XFORM-011): counterclockwise positive. */
    readonly degrees: '{n}°';
    /** The dialogs (UI-018, A11Y-011). */
    readonly close: 'Close';
    readonly cancel: 'Cancel';
    readonly apply: 'Apply';
    readonly add: 'Add';
    readonly atomProperties: 'Atom Properties';
    /** Atom Properties' fields (ATOM-029). */
    readonly fieldElement: 'Element';
    readonly fieldNumber: 'Atomic number';
    readonly fieldAlias: 'Alias';
    readonly fieldCharge: 'Charge';
    readonly fieldIsotope: 'Isotope';
    readonly fieldValence: 'Valence';
    readonly fieldRadical: 'Radical';
    /** The radical options, as Ketcher names them, never as valences (fix round 4). */
    readonly radicalNone: 'None';
    readonly radicalMono: 'Monoradical';
    readonly radicalDi: 'Diradical';
    /** More radical electrons than the user may choose, shown only as an atom has them (`[N]`, `[C]`). */
    readonly radicalTri: 'Triradical';
    readonly radicalTetra: 'Tetraradical';
    readonly radicalElectrons: {
        readonly one: '{count} radical electron';
        readonly other: '{count} radical electrons';
    };
    /** Its messages (ATOM-031, UI-018). */
    readonly notAnElement: '"{typed}" is not an element symbol';
    readonly elementNeeded: 'Type an element symbol';
    readonly chargeRange: 'A whole number from -15 to 15';
    readonly isotopeRange: 'A whole number from 1 to 999, or blank';
    readonly valenceRange: 'A whole number from 0 to 8, or blank';
    /** Beside the alias field, and what refuses a longer label or alias from a script (fix round 3; ALIAS_LIMIT). */
    readonly aliasLimit: 'An alias is at most 80 characters';
    /** A field's hint while what it holds is not valid, shown as the error it is (fix round 4; WCAG 1.4.1). */
    readonly notValid: 'Not valid: {hint}';
    /** The periodic table (ATOM-009). */
    readonly periodicElements: 'Elements';
    readonly detailsNumber: 'Atomic number {n}';
    readonly detailsMass: 'Mass {mass}';
    readonly nothingChosen: 'Choose an element';
    /**
     * The context menus (spike context-menus, C2, C10): each named after what it acts on; their items, dialog titles with an
     * ellipsis where they open a dialog, the rest in sentence case; the groups of radio items. Bond types are named as the
     * toolbar names them (singleBond …), and Copy As, Paste Special, Paste, Cut, Copy, Clean Up and Fit to view as theirs.
     */
    readonly menuAtom: 'Atom';
    readonly menuBond: 'Bond';
    readonly menuSelection: 'Selection';
    readonly menuCanvas: 'Canvas';
    readonly editLabel: 'Edit label';
    readonly opensDialog: '{title}…';
    readonly selectStructure: 'Select structure';
    readonly selectAll: 'Select all';
    readonly changeDirection: 'Change direction';
    /** The atom's and the selection's Charge submenu (spike charges-checks, K4; context-menus' C3): the charge tools' names. */
    readonly charge: 'Charge';
    /** A charge the menu would take past ±15 (CHG-006), as the status line says it: "No charge past +15 on this N". */
    readonly chargeLimitUp: 'No charge past +15 on this {symbol}';
    readonly chargeLimitDown: 'No charge past −15 on this {symbol}';
    readonly deleteItem: 'Delete';
    readonly bondType: 'Bond type';
    readonly doubleBondPosition: 'Double bond position';
    readonly positionLeft: 'Left';
    readonly positionCentre: 'Centre';
    readonly positionRight: 'Right';
    /** Bond Properties (spike context-menus, C8; BOND-018): its title and its one field. */
    readonly bondProperties: 'Bond Properties';
    readonly fieldType: 'Type';
    /** Atom Properties over several atoms (C9, ATOM-030): what a field shows where the atoms differ. */
    readonly severalValues: 'Several values';
    /**
     * Spike stereo (docs/decisions/stereo.md, S5 to S10): a double bond's E/Z in its menu; the canvas menu's check boxes;
     * the Enhanced Stereochemistry dialog (Ketcher's and ChemDraw's title), its types and groups; the flag, in ChemDraw's
     * words (P4); a 180° rotation refused across a double bond with E or Z (S6, P9). The groups' names, &1 and or1, are
     * notations, never translated.
     */
    readonly ezGroup: 'E/Z';
    readonly ezE: 'E';
    readonly ezZ: 'Z';
    readonly ezUnknown: 'Unknown';
    readonly chiralFlag: 'Chiral flag';
    readonly showStereochemistry: 'Show stereochemistry';
    readonly enhancedStereo: 'Enhanced Stereochemistry';
    /**
     * Milestone-3-fixes (F11, the product owner's call of 2026-10-03): a potential stereocentre's Stereochemistry submenu,
     * R, S or Unspecified (CIP labels, as E and Z are), and why Enhanced Stereochemistry… is greyed out.
     */
    readonly stereochemistry: 'Stereochemistry';
    readonly stereoR: 'R';
    readonly stereoS: 'S';
    readonly stereoUnspecifiedItem: 'Unspecified';
    readonly setRorSFirst: 'Set R or S first';
    /**
     * Spike abbreviations (docs/decisions/abbreviations.md, A5 to A8; I18N-006: the abbreviations' names are notations, never
     * translated, and are passed into these as they are): a contracted abbreviation's menu, named "Abbreviation", and its
     * items, ChemDraw's and Ketcher's names (title case, as theirs); the prompt before an edit inside an expanded abbreviation
     * (P5, Ketcher's), its title, question and buttons; the group keys' row of the shortcuts sheet (KEY-009); what an undo of
     * Contract Label says.
     */
    readonly menuAbbreviation: 'Abbreviation';
    readonly expandAbbreviation: 'Expand Abbreviation';
    readonly contractAbbreviation: 'Contract Abbreviation';
    readonly removeAbbreviation: 'Remove Abbreviation';
    readonly contractLabel: 'Contract Label';
    readonly dissolveTitle: 'Remove abbreviation';
    readonly dissolveQuestion: {
        readonly one: 'Remove abbreviation {labels}?';
        readonly other: 'Remove abbreviations {labels}?';
    };
    readonly dissolveHint: 'This edit changes atoms of the abbreviation. Its atoms and bonds stay, drawn in full.';
    readonly dissolveRemove: 'Remove and edit';
    readonly sheetAbbreviation: 'Abbreviation {label}';
    /** A centre the engine cannot give that configuration (a cage's): "No R configuration for this C". */
    readonly noCentreConfig: 'No {config} configuration for this {symbol}';
    readonly stereoType: 'Group';
    readonly stereoAbs: 'Absolute (abs)';
    readonly stereoAnd: 'AND';
    readonly stereoOr: 'OR';
    readonly stereoNone: 'None';
    readonly andGroup: 'AND group';
    readonly orGroup: 'OR group';
    readonly newGroup: 'New ({name})';
    readonly stereoTargets: {
        readonly one: '{count} stereocentre';
        readonly other: '{count} stereocentres';
    };
    readonly flagAbsolute: 'Absolute';
    readonly flagRacemic: 'Racemic';
    readonly flagRelative: 'Relative';
    readonly cannotTurnEz: 'Cannot rotate by 180° across a double bond with E or Z: it would change its E/Z';
    /**
     * The S-Group Properties dialog (spike sgroups, G4; Ketcher's title, fields and words), its checks and the engine's refusals
     * (ABBR-033, 035, 037, 040: "Partial S-group overlapping is not allowed." is Ketcher's), the menu's items (G5), what an undo
     * names, and a group's tooltip (ABBR-043). Labels, field names, values and types are passed in as they are.
     */
    readonly sgroupDialog: 'S-Group Properties';
    readonly sgroupType: 'Type';
    readonly sgroupTypeSup: 'Superatom';
    readonly sgroupTypeDat: 'Data';
    readonly sgroupTypeMul: 'Multiple group';
    readonly sgroupTypeSru: 'SRU polymer';
    readonly sgroupName: 'Name';
    readonly sgroupContext: 'Context';
    readonly contextAtom: 'Atom';
    readonly contextBond: 'Bond';
    readonly contextFragment: 'Fragment';
    readonly contextGroup: 'Group';
    readonly contextMultifragment: 'Multifragment';
    readonly sgroupFieldName: 'Field name';
    readonly sgroupFieldValue: 'Field value';
    readonly sgroupPlacement: 'Placement';
    readonly placementAbsolute: 'Absolute';
    readonly placementRelative: 'Relative';
    readonly placementAttached: 'Attached';
    readonly sgroupCount: 'Repeat count';
    readonly sgroupLabel: 'Polymer label';
    readonly sgroupConnect: 'Repeat pattern';
    readonly connectHT: 'Head-to-tail';
    readonly connectHH: 'Head-to-head';
    readonly connectEU: 'Either unknown';
    readonly sgroupNeedsName: 'Enter a name';
    readonly sgroupNameTooLong: 'A name is at most {n} characters';
    readonly sgroupNeedsField: 'Enter a field name and a value';
    readonly sgroupFieldNameTooLong: 'A field name is at most {n} characters';
    readonly sgroupValueTooLong: 'A value is at most {n} characters a line';
    readonly sgroupCountRange: 'The repeat count is a whole number from 1 to {n}';
    readonly sgroupNeedsLabel: 'Enter a polymer label, without quotation marks';
    readonly sgroupLabelTooLong: 'A polymer label is at most {n} characters';
    readonly sgroupPartialOverlap: 'Partial S-group overlapping is not allowed.';
    readonly sgroupContextMismatch: 'What was chosen does not fit the context {context}';
    readonly sgroupMulBonds: 'A multiple group needs no crossing bond, or two';
    readonly sgroupRefused: 'This S-group cannot be made here';
    readonly editSGroup: 'Edit S-Group';
    readonly removeSGroup: 'Remove S-Group';
    readonly stepSGroup: 'S-group';
    readonly sgroupTipData: 'Data S-group: {name} = {value}';
    readonly sgroupTipSru: 'SRU polymer: {label} ({connect})';
    readonly sgroupTipMul: 'Multiple group: {count}';
    readonly sgroupTipSup: 'Superatom: {label}';
    readonly sgroupTipOther: '{type} S-group';
    /** The R-Group dialog (RGROUP-001, RGROUP-003; Ketcher's title). Its toggles read R1 to R32, R labels never translated. */
    readonly ok: 'OK';
    readonly rgroupDialog: 'R-Group';
    readonly rgroupNumbers: 'R-group numbers';
    readonly rgroupHint: 'Choose one or more; none makes an R label a carbon again';
    /** The Attachment Points dialog (RGROUP-014; Ketcher's title and check boxes). */
    readonly attachmentDialog: 'Attachment Points';
    readonly primaryPoint: 'Primary attachment point';
    readonly secondaryPoint: 'Secondary attachment point';
    /**
     * What screen readers hear (spike keyboard-a11y, K3; A11Y-007, A11Y-009, I18N-004): counts, the canvas's description,
     * each edit, each move of the keyboard's hotspot, each selection, a tool or a zoom chosen by key. An atom is named by
     * its label as drawn (a symbol) and its number counted from 1; a bond by its number and its type.
     */
    readonly atomCount: {
        readonly one: '{count} atom';
        readonly other: '{count} atoms';
    };
    readonly bondCount: {
        readonly one: '{count} bond';
        readonly other: '{count} bonds';
    };
    readonly both: '{first} and {second}';
    /** The canvas's description (A11Y-007): what is drawn, kept current. The formula joins it with spike 7. */
    readonly descriptionEmpty: 'Empty';
    readonly description: '{counts}. SMILES: {smiles}';
    readonly withFormula: '{formula}, {counts}';
    /** The keyboard's hotspot moved (K1): onto an atom, or a bond; or nothing lies that way, and it stayed. */
    readonly hotAtom: '{label}, atom {number}, {bonds}';
    readonly hotBond: 'Bond {number}, {type}, between atoms {first} and {second}';
    readonly nothingThatWay: 'Nothing that way';
    /** A bond's type as an announcement says it. */
    readonly typeSingle: 'single';
    readonly typeDouble: 'double';
    readonly typeTriple: 'triple';
    readonly typeAromatic: 'aromatic';
    readonly typeWedge: 'wedge';
    readonly typeHash: 'hashed wedge';
    readonly typeWavy: 'wavy';
    /** A crossed double bond (milestone-3-fixes, F4), as "Bond 2 crossed" says it. */
    readonly typeCrossed: 'crossed';
    readonly typeOther: 'other';
    /** The selection changed, by any input (K3). */
    readonly selected: '{what} selected';
    readonly selectionCleared: 'Selection cleared';
    /** An edit, one message per change event (K3). */
    readonly added: 'Added {what}';
    readonly deleted: 'Deleted {what}';
    readonly pasted: 'Pasted {what}';
    readonly canvasCleared: 'Canvas cleared';
    readonly bondNow: 'Bond {number} {type}';
    readonly bondsNow: '{bonds} {type}';
    readonly atomNow: 'Atom {number} {label}';
    readonly atomsNow: '{atoms} {label}';
    readonly changed: 'Changed {what}';
    readonly moved: 'Moved {what}';
    readonly rotated: 'Rotated {what}';
    readonly flipped: 'Flipped {what}';
    readonly joined: 'Joined {what}';
    readonly laidOut: 'Laid out the drawing';
    readonly cleanedUp: 'Cleaned up the drawing';
    /** Milestone-3-fixes (F4, F3, F11): a Clean Up of a selection, Straighten, a check box, a centre set unspecified, a tool chosen to draw. */
    readonly cleanedAtoms: 'Cleaned up {what}';
    readonly straightened: 'Straightened';
    readonly switchedOn: '{name} on';
    readonly switchedOff: '{name} off';
    readonly stereoUnspecified: 'unspecified';
    readonly toolThen: '{tool}. {what}';
    readonly undone: 'Undo: {name}';
    readonly redone: 'Redo: {name}';
    readonly newDrawing: 'New drawing: {what}';
    /** The names of the edits, as an undo or a redo says them ("Undo: Paste"). */
    readonly stepAdd: 'Add';
    readonly stepDelete: 'Delete';
    readonly stepBond: 'Bond change';
    readonly stepLabel: 'Label';
    readonly stepRing: 'Ring';
    readonly stepTemplate: 'Template';
    readonly stepMove: 'Move';
    readonly stepDuplicate: 'Duplicate';
    readonly stepJoin: 'Join';
    readonly stepRotate: 'Rotate';
    readonly stepFlip: 'Flip';
    readonly stepEdit: 'Edit';
    /** A tool or a zoom chosen by key (K3). */
    readonly toolChosen: '{tool} tool';
    readonly zoomTo: 'Zoom {percent}';
    /** A notice's close button (UI-020): it stays until the next edit, Escape, this, or 6 s. */
    readonly dismiss: 'Dismiss';
    /**
     * Spike settings-touch (K6): the top bar's last group and its gear; the Settings dialog (Ketcher's title), a row for each
     * of the sketcher's options, its values, and Reset.
     */
    readonly groupSketcher: 'Sketcher';
    readonly settings: 'Settings';
    /** The gear as one menu button (spike bar-fit, B9): its name and its menu's, and its items, which open dialogs. */
    readonly optionsMenu: 'Settings and help';
    readonly settingsItem: 'Settings…';
    readonly keysItem: 'Keyboard shortcuts…';
    /** The gear's check item that switches query mode (spike reports-oct5, R3), and what is said when a user switches it. */
    readonly queryModeItem: 'Query mode';
    readonly queryModeOn: 'Query mode on';
    readonly queryModeOff: 'Query mode off';
    readonly settingAtomColours: 'Colour atoms by element';
    readonly settingWarningMarks: 'Mark atoms with problems';
    readonly settingReadout: 'Formula and mass';
    readonly readoutShown: 'Shown';
    readonly readoutCollapsed: 'Collapsed';
    readonly readoutHidden: 'Hidden';
    readonly settingStereoLabels: 'Show R, S, E and Z labels';
    readonly settingStereoFlag: 'Show the stereo flag';
    readonly settingUnspecified: 'Draw unspecified double bonds crossed';
    readonly settingMethyl: 'Label terminal carbons CH₃';
    readonly settingKeymap: 'Keys';
    readonly keymapChemDraw: 'ChemDraw';
    readonly keymapKetcher: 'Ketcher';
    readonly reset: 'Reset';
    /**
     * The shortcuts sheet (K7, P13; KEY-006): its button and title, its sections by where the keys act, its columns, and the
     * names of what keys do where the toolbars have none.
     */
    readonly keyboardShortcuts: 'Keyboard shortcuts';
    readonly sheetWhat: 'Action';
    readonly sheetKeys: 'Keys';
    readonly sheetAnywhere: 'Anywhere in the sketcher';
    readonly sheetNothingHovered: 'With nothing under the pointer';
    readonly sheetOverAtom: 'Over an atom';
    readonly sheetOverBond: 'Over a bond';
    readonly sheetSelection: 'With a selection';
    readonly sheetInTurn: '{tools}, in turn';
    readonly sheetFuse: 'Fuse {ring}';
    readonly sheetNudge: 'Move the selection {n} pt';
    readonly sheetRotate: 'Rotate the selection {n}°';
    readonly sheetCopyAs: 'Copy As {format}';
    readonly sheetPasteAs: 'Paste Special {format}';
    readonly sheetPositionLeft: 'Double bond on the left';
    readonly sheetPositionCentre: 'Double bond centred';
    readonly sheetPositionRight: 'Double bond on the right';
    readonly sheetDeuterium: 'Deuterium';
    readonly sheetUnlabel: 'Carbon, or delete a carbon';
    readonly sheetSelectNone: 'Deselect all';
    readonly sheetSelectLast: 'Select the structure drawn last';
    readonly sheetTurnHorizontal: 'Rotate 180° horizontally';
    readonly sheetTurnVertical: 'Rotate 180° vertically';
    readonly sheetMenu: 'Context menu';
    readonly sheetWalk: 'Move the keyboard hotspot';
    readonly sheetJump: 'Move the hotspot three bonds';
    readonly sheetHop: 'Move the hotspot to the next structure';
    readonly sheetStart: 'Draw, or put the hotspot on the drawing';
    readonly sheetEndSelection: 'End the selection';
    /** What the number keys and a sprout over an atom (KEY-011), as the sheet names them. */
    readonly sprout1: 'A bond';
    readonly sprout2: 'C=O, or an acetyl';
    readonly sprout3: 'A phenyl';
    readonly sprout4: 'A wedged bond';
    readonly sprout5: 'A hashed wedge bond';
    readonly sprout6: 'A cyclohexane';
    readonly sprout7: 'A cyclopentane';
    readonly sprout8: '=CH₂';
    readonly sprout9: 'Two methyls, or an isopropyl';
    readonly sprout0: 'A long bond';
    /** Small boxes (K8): a bar's scroll buttons. */
    readonly scrollBack: 'Scroll back';
    readonly scrollForward: 'Scroll forward';
    /** The query bond tools (QUERY-006), Ketcher's names, and their palette. */
    readonly anyBond: 'Any bond';
    readonly singleOrDoubleBond: 'Single or double bond';
    readonly singleOrAromaticBond: 'Single or aromatic bond';
    readonly doubleOrAromaticBond: 'Double or aromatic bond';
    readonly queryBonds: 'Query bonds';
    /** The query atom tool (Q8): its control's name, and its button's for what it holds. */
    readonly queryAtom: 'Query atom';
    readonly queryAtomList: 'Atom list {label}';
    readonly queryAtomNotList: 'NOT list {label}';
    readonly queryAtomGeneric: 'Generic atom {label}';
    /** Query properties ▸ (Q5) and Atom Properties' Query section (Q6): the fields' names, then their choices. */
    readonly queryProperties: 'Query properties';
    readonly queryRingBondCount: 'Ring bond count';
    readonly queryHCount: 'H count';
    readonly querySubstitution: 'Substitution count';
    readonly queryUnsaturated: 'Unsaturated';
    readonly queryAromaticity: 'Aromaticity';
    readonly queryImplicitH: 'Implicit H count';
    readonly queryRingMembership: 'Ring membership';
    readonly queryRingSize: 'Ring size';
    readonly queryConnectivity: 'Connectivity';
    readonly queryAny: 'Any';
    readonly queryAsDrawn: 'As drawn';
    readonly queryNone: 'None';
    readonly queryOrMore: '{n} or more';
    readonly queryAromatic: 'Aromatic';
    readonly queryAliphatic: 'Aliphatic';
    readonly clearQueryProperties: 'Clear query properties';
    /** Why the query items are greyed out while a custom query is set (QUERY-015). */
    readonly customAtomQuerySet: 'A custom query is set: change it in Atom Properties';
    readonly customBondQuerySet: 'A custom query is set: change it in Bond Properties';
    /** Topology ▸ and Bond Properties' Topology (QUERY-007). */
    readonly topology: 'Topology';
    readonly topologyEither: 'Either';
    readonly topologyRing: 'Ring';
    readonly topologyChain: 'Chain';
    /** Atom Properties' and Bond Properties' Query sections (Q6, Q7). */
    readonly querySection: 'Query';
    readonly customQuery: 'Custom query';
    readonly customQuerySmarts: 'Custom query SMARTS';
    readonly smartsNotRead: 'Not valid SMARTS at character {position}: {message}';
    readonly listTooLong: 'An atom list holds at most 16 elements';
    readonly listEmpty: 'An atom list needs at least one element';
    /** The periodic table in query mode (Q8, QUERY-004): its modes and its Generics row, RDKit's meanings. */
    readonly tableMode: 'Pick';
    readonly tableSingle: 'Single';
    readonly tableList: 'List';
    readonly tableNotList: 'Not list';
    readonly genericsRow: 'Generics';
    readonly genericA: 'A: any atom but hydrogen';
    readonly genericAH: 'AH: any atom';
    readonly genericQ: 'Q: any atom but carbon and hydrogen';
    readonly genericQH: 'QH: any atom but carbon';
    readonly genericX: 'X: a halogen';
    readonly genericXH: 'XH: a halogen or hydrogen';
    readonly genericM: 'M: a metal';
    readonly genericMH: 'MH: a metal or hydrogen';
    readonly genericStar: '*: any atom';
    /** SMARTS read in query mode only (E5, QUERY-016; P1). */
    readonly queryModeOnly: '{name}: in query mode only';
    readonly pasteQueryOnly: 'Pasting {format} needs query mode';
    readonly dropQueryOnly: 'Dropping {format} files needs query mode';
    /** A query element the host turned off (API-028). */
    readonly queryElementOff: '{name}: turned off in this sketcher';
    /** The shortcuts sheet's row for x over an atom (KEY-008, P7). */
    readonly sheetGenericX: 'The generic atom X (query mode)';
};
/** A message's name. */
export type StringKey = keyof typeof STRINGS;
/** The messages with placeholders or plural forms: read through the catalog's functions, never as text. */
type FormattedKey = 'groupCount' | 'percent' | 'zoomButton' | 'withHotkey' | 'notAvailableYet' | 'noRoomBond' | 'noRoomRing' | 'noRoomChain' | 'noRoomTemplate' | 'noRoomMerge' | 'thisAtom' | 'noRoomForRadical' | 'problemValence' | 'chargeLimitUp' | 'chargeLimitDown' | 'tooClose' | 'cannotFlip' | 'cannotTurn' | 'noFreeAttachment' | 'pasteNotYet' | 'clipboardBlocked' | 'copyAsBlocked' | 'copiedLost' | 'copiedV3000' | 'copiedV3000Lost' | 'copiedV3000Needed' | 'copiedNone' | 'lossAlias' | 'lossRLabelNumbers' | 'lossAttachmentPoints' | 'lossStereoGroups' | 'lossChiralFlag' | 'lossQuery' | 'lossSgroup' | 'lossHydrogenBond' | 'lossZeroOrder' | 'lossAtomMaps' | 'lossStereo' | 'lossRLabel' | 'lossDummy' | 'lossDative' | 'lossOther' | 'pasteNotAs' | 'pasteNotAsMol' | 'pasteSpecialHint' | 'notReadAs' | 'notReadAsMol' | 'dropNothing' | 'dropNotYet' | 'pasteNoText' | 'readPasted' | 'readDropped' | 'readWritten' | 'readLeftOut' | 'readAlso' | 'leftStructure' | 'leftText' | 'leftGraphic' | 'leftArrow' | 'leftPicture' | 'leftBracket' | 'keptAliases' | 'lostRecords' | 'lostSdData' | 'lostRgroupDefinitions' | 'lostCoordinates3d' | 'lostSgroup' | 'lostQuery' | 'lostCxSection' | 'lostMolfile' | 'lostOther' | 'degrees' | 'radicalElectrons' | 'notAnElement' | 'notValid' | 'detailsNumber' | 'detailsMass' | 'opensDialog' | 'atomCount' | 'bondCount' | 'both' | 'description' | 'withFormula' | 'hotAtom' | 'hotBond' | 'selected' | 'added' | 'deleted' | 'pasted' | 'bondNow' | 'bondsNow' | 'atomNow' | 'atomsNow' | 'changed' | 'moved' | 'rotated' | 'flipped' | 'joined' | 'undone' | 'redone' | 'newDrawing' | 'noCentreConfig' | 'cleanedAtoms' | 'switchedOn' | 'switchedOff' | 'toolThen' | 'toolChosen' | 'newGroup' | 'stereoTargets' | 'zoomTo' | 'sheetInTurn' | 'sheetFuse' | 'sheetNudge' | 'sheetRotate' | 'sheetCopyAs' | 'sheetPasteAs' | 'dissolveQuestion' | 'sheetAbbreviation' | 'sgroupNameTooLong' | 'sgroupFieldNameTooLong' | 'sgroupValueTooLong' | 'sgroupCountRange' | 'sgroupLabelTooLong' | 'sgroupContextMismatch' | 'sgroupTipData' | 'sgroupTipSru' | 'sgroupTipMul' | 'sgroupTipSup' | 'sgroupTipOther' | 'queryAtomList' | 'queryAtomNotList' | 'queryAtomGeneric' | 'queryOrMore' | 'smartsNotRead' | 'queryModeOnly' | 'pasteQueryOnly' | 'dropQueryOnly' | 'queryElementOff';
/** The messages read as plain text. */
export type PlainKey = Exclude<StringKey, FormattedKey>;
/** What a refused gesture was for, as a refusal says it (F3). */
export type RefusalWhy = 'bond' | 'ring' | 'chain' | 'template' | 'merge';
/** Parameters of a template: each `{name}` is replaced by its value, written as it is. */
export type Params = Readonly<Record<string, string | number>>;
/**
 * A sketcher's strings, in its locale (K5): the catalog's plain messages as text, and the rest as functions. Each
 * sketcher has its own, so that two on one page can speak two languages (API-032).
 */
export type Catalog = {
    readonly [K in PlainKey]: string;
} & {
    /** The locale in force, as `Intl` writes it (`en`, `de-CH`, `en-XA`). */
    readonly locale: string;
    /** A message by name, its placeholders filled. */
    format(key: StringKey, params?: Params): string;
    /** A plural message for `count`, in the form the locale takes for it, `{count}` and the other placeholders filled. */
    count(key: StringKey, count: number, params?: Params): string;
    groupCount(group: string, count: number): string;
    percent(n: number): string;
    withHotkey(name: string, key: string): string;
    /** The tooltip around its key, which the toolbars set in a `<kbd>` (fix round 4): the text before the key, and after it. */
    aroundHotkey(name: string): readonly [string, string];
    notAvailableYet(name: string): string;
    /** An element's name (ATOM-001): the host's `element.<Symbol>`, or the element data's English name. */
    elementName(symbol: string): string;
    noRoom(atom: string, what: RefusalWhy): string;
    /** An atom on the canvas, as a refusal names it: "this C" (K14). */
    thisAtom(symbol: string): string;
    noRoomForRadical(symbol: string): string;
    problemValence(symbol: string): string;
    chargeLimit(symbol: string, delta: 1 | -1): string;
    tooClose(atom: string): string;
    cannotFlip(bonds: number, how: 'mirror' | 'turn'): string;
    noFreeAttachment(atom: string): string;
    pasteNotYet(format: string): string;
    clipboardBlocked(action: string, key: string): string;
    /** A format's name, a notation, never translated. */
    formatName(format: string): string;
    copyAsBlocked(format: string): string;
    /** A Copy As's losses (milestone-4-fixes, F4), each named by `lossName`. */
    copiedLost(format: string, what: readonly string[]): string;
    /** MOL V2000 asked for and V3000 written (F4): why, and what V3000 did not write. */
    copiedV3000(why: readonly string[], what: readonly string[]): string;
    /** A loss table code's words (F4): "enhanced stereo groups"; any other code, "some features". */
    lossName(code: string): string;
    /** Copy As of a part `format` cannot hold at all (spike inchi, N6): "Not copied: InChI cannot hold R labels". */
    copiedNone(format: string, what: readonly string[]): string;
    pasteNotAs(format: string): string;
    pasteSpecialHint(key: string): string;
    notReadAs(format: string): string;
    dropNothing(name: string): string;
    dropNotYet(format: string): string;
    /** A paste with no text (spike cdxml, X5), naming ChemDraw's Copy As CDXML key as this sketcher writes it. */
    pasteNoText(key: string): string;
    /** What a CDXML read left out of one kind (X3): "1 text", "2 arrows", "1 structure that could not be read". */
    leftOut(kind: 'structure' | 'text' | 'graphic' | 'arrow' | 'picture' | 'bracket', n: number): string;
    /** The read report (X3): "Pasted 2 structures; left out: 1 text, 1 arrow; 1 label kept as an alias". */
    readReport(how: 'paste' | 'drop' | 'write', structures: number, leftOut: readonly string[], aliases: number): string;
    /** What a read dropped, by the engine's code (milestone-4-fixes, F1): "the records after the first"; any other code, "some of the text". */
    readLost(code: string): string;
    degrees(n: number): string;
    radicalMore(n: number): string;
    notAnElement(typed: string): string;
    notValid(hint: string): string;
    detailsNumber(n: number): string;
    detailsMass(mass: number): string;
    opensDialog(title: string): string;
    /** R1 … R32: R labels, never translated. */
    rName(n: number): string;
    /** Enhanced stereo groups' names, &1 and or1: notations, never translated (spike stereo). */
    andName(n: number): string;
    orName(n: number): string;
    newGroup(name: string): string;
    stereoTargets(n: number): string;
    /** The prompt before an edit inside expanded abbreviations (spike abbreviations, P5): "Remove abbreviation OMe?". */
    dissolveQuestion(labels: readonly string[]): string;
    /** The S-Group Properties dialog's checks and refusals (spike sgroups, G4), with their limits and the context's name. */
    sgroupNameTooLong(n: number): string;
    sgroupFieldNameTooLong(n: number): string;
    sgroupValueTooLong(n: number): string;
    sgroupCountRange(n: number): string;
    sgroupLabelTooLong(n: number): string;
    sgroupContextMismatch(context: string): string;
    /** A type's, a context's, a placement's and a repeat pattern's names in the dialog: "Multiple group", "Fragment", "Head-to-tail". */
    sgroupTypeName(type: string): string;
    sgroupContextName(context: string): string;
    sgroupPlacementName(placement: string): string;
    sgroupConnectName(connect: string): string;
    /** A group's tooltip (ABBR-043): "Data S-group: pKa = 4.2", "SRU polymer: n (Head-to-tail)", "Multiple group: 3". */
    sgroupTipData(name: string, value: string): string;
    sgroupTipSru(label: string, connect: string): string;
    sgroupTipMul(count: string): string;
    sgroupTipSup(label: string): string;
    sgroupTipOther(type: string): string;
    /** "3 atoms", "1 bond". */
    atoms(n: number): string;
    bonds(n: number): string;
    /** "3 atoms and 2 bonds", leaving out a part that is none (but both none: "0 atoms"). */
    atomsAndBonds(atoms: number, bonds: number): string;
};
/** A template with its placeholders filled; a placeholder with no value is left as it is. */
export declare function fill(template: string, params?: Params): string;
/** The pseudo-locale, built in: English accented, about 40% longer, bracketed (K5). */
export declare const PSEUDO_LOCALE = "en-XA";
/**
 * English text as the pseudo-locale writes it (I18N-003): its letters accented, `{placeholders}` kept as they are, about
 * 40% longer with tildes, in brackets: "Single bond" is "[Šîñĝļé ƀöñð~~~~]". Text with no letter outside its
 * placeholders ("{n}%", "{name} ({key})") is left as it is: it says nothing a translation would change.
 */
export declare function pseudo(text: string): string;
/**
 * The locale a host asked for, as `Intl` writes it (`de-ch` is `de-CH`); English for anything `Intl` cannot read, and for
 * nothing given. The pseudo-locale is its own.
 */
export declare function toLocale(value: unknown): string;
/** A host's `strings`, kept to what the catalog can take: known names (and element names) with messages. */
export declare function toStrings(value: unknown): SketcherStrings;
/**
 * A sketcher's catalog (K5): English, overlaid by the pseudo-locale when it is `en-XA`, overlaid by the host's `strings`
 * (I18N-002). Plurals follow the locale's `Intl.PluralRules` (I18N-004).
 */
export declare function catalogOf(localeAsked?: unknown, stringsGiven?: unknown): Catalog;
/** English, the default: every sketcher's catalog until its host gives it a locale or strings. */
export declare const EN: Catalog;
export {};
