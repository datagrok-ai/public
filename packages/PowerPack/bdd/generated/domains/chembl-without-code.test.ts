/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/domains/chembl-without-code.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {childRows, columnSemType, followLink, goBack, leaveOut, openRow, pickTable, schemaGone, settingOn, tableCount} from '../../bindings/domains.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, isExpanded, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, connectionsOnServer, dialogCloses, ownConnectionGone, shellSettingPutBack, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Connect ChEMBL and get an application over it without code", () => {
  const session = feature(test, "features/domains/chembl-without-code.feature", import.meta.url);
  test("Connect, bind every table, and walk ChEMBL along its own links", {tag: ["@guide", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And simple mode is off", () => simpleModeOff(page));
    await session.step(12, "And the user's own connection \"ChEMBL37\" is removed now and at feature end", () => ownConnectionGone(page, "ChEMBL37"));
    await session.step(13, "And the domain schema \"chembl37\" is removed now and at feature end", () => schemaGone(page, "chembl37"));
    await session.step(14, "And the \"enableDomainDatabases\" shell setting is put back at feature end", () => shellSettingPutBack(page, "enableDomainDatabases"));
    await session.step(15, "And the \"enableDomainDatabases\" shell setting is on", () => settingOn(page, "enableDomainDatabases"));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(19, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")), undefined, "A Postgres connection to ChEMBL, with a read-only login");
    await session.step(20, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(21, "When user enters \"ChEMBL37\" into Name input in \"Add new connection\" dialog", () => enterInto(page, "ChEMBL37", el("Name input in \"Add new connection\" dialog")));
    await session.step(22, "And user enters \"db.datagrok.ai\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "db.datagrok.ai", el("Server input in \"Add new connection\" dialog")));
    await session.step(23, "And user enters \"54325\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "54325", el("Port input in \"Add new connection\" dialog")));
    await session.step(24, "And user enters \"chembl\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "chembl", el("Db input in \"Add new connection\" dialog")));
    await session.step(25, "And user enters \"datagrok\" into Login input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("Login input in \"Add new connection\" dialog")));
    await session.step(26, "And user enters \"datagrok\" into Password input in \"Add new connection\" dialog", () => enterInto(page, "datagrok", el("Password input in \"Add new connection\" dialog")));
    await session.step(27, "And user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(28, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(29, "And 1 connection named \"ChEMBL37\" should be on the server", () => connectionsOnServer(page, 1, "ChEMBL37"));
    await session.step(31, "When user expands Databases---Postgres tree node inside browse tree", () => expand(page, el("Databases---Postgres tree node inside browse tree")), undefined, "Right-click its schema › Create Domain Schema…");
    await session.step(32, "And user expands Databases---Postgres---ChEMBL37 tree node inside browse tree", () => expand(page, el("Databases---Postgres---ChEMBL37 tree node inside browse tree")));
    await session.step(33, "And user expands Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree", () => expand(page, el("Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree")));
    await session.step(34, "And user picks \"Create Domain Schema...\" from the context menu of Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree", () => pickFromContextMenu(page, "Create Domain Schema...", el("Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree")));
    await session.step(35, "Then \"Create domain schema\" dialog should be visible", () => shouldBe(page, el("\"Create domain schema\" dialog"), "visible"));
    await session.step(37, "And \"Create domain schema\" dialog should contain text \"71 of 71 bindable tables\"", () => shouldContainText(page, el("\"Create domain schema\" dialog"), "71 of 71 bindable tables"), undefined, "All 71 tables are bindable: keys and links come from the database, structures are recognized");
    await session.step(39, "When user enters \"ChEMBL\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "ChEMBL", el("\"Friendly name\" input in \"Create domain schema\" dialog")), undefined, "Name the schema, and give the tables people will open a readable name and a caption column");
    await session.step(40, "And user enters \"chembl37\" into Identifier input in \"Create domain schema\" dialog", () => enterInto(page, "chembl37", el("Identifier input in \"Create domain schema\" dialog")));
    await session.step(41, "And user leaves the \"chembl_id_lookup\" table out of schema design", () => leaveOut(page, "chembl_id_lookup"));
    await session.step(42, "And user clicks on the \"target_dictionary\" table in schema design", () => pickTable(page, "target_dictionary"));
    await session.step(43, "And user enters \"Targets\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Targets", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(44, "And user selects \"pref_name\" in \"Name column\" input in \"Create domain schema\" dialog", () => selectIn(page, "pref_name", el("\"Name column\" input in \"Create domain schema\" dialog")));
    await session.step(45, "And user clicks on the \"molecule_dictionary\" table in schema design", () => pickTable(page, "molecule_dictionary"));
    await session.step(46, "And user enters \"Compounds\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Compounds", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(47, "And user selects \"chembl_id\" in \"Name column\" input in \"Create domain schema\" dialog", () => selectIn(page, "chembl_id", el("\"Name column\" input in \"Create domain schema\" dialog")));
    await session.step(48, "And user selects \"pref_name\" in \"Searchable\" input in \"Create domain schema\" dialog", () => selectIn(page, "pref_name", el("\"Searchable\" input in \"Create domain schema\" dialog")));
    await session.step(49, "And user clicks on the \"assays\" table in schema design", () => pickTable(page, "assays"));
    await session.step(50, "And user enters \"Assays\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Assays", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(51, "And user selects \"chembl_id\" in \"Name column\" input in \"Create domain schema\" dialog", () => selectIn(page, "chembl_id", el("\"Name column\" input in \"Create domain schema\" dialog")));
    await session.step(52, "And user clicks on the \"docs\" table in schema design", () => pickTable(page, "docs"));
    await session.step(53, "And user enters \"Documents\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Documents", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(54, "And user selects \"chembl_id\" in \"Name column\" input in \"Create domain schema\" dialog", () => selectIn(page, "chembl_id", el("\"Name column\" input in \"Create domain schema\" dialog")));
    await session.step(55, "And user clicks on \"Next\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Next\" button in \"Create domain schema\" dialog")));
    await session.step(57, "And user clicks on \"Validate\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Validate\" button in \"Create domain schema\" dialog")), undefined, "Validate checks every table against the live database before anything is registered");
    await session.step(58, "Then \"Create domain schema\" dialog should contain text \"No findings\"", () => shouldContainText(page, el("\"Create domain schema\" dialog"), "No findings"));
    await session.step(59, "When user clicks on \"Create\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Create\" button in \"Create domain schema\" dialog")));
    await session.step(60, "Then the \"Create domain schema\" dialog should close", () => dialogCloses(page, "Create domain schema"));
    await session.step(61, "And the \"chembl37.target_dictionary\" domain table should hold 18552 rows on the server", () => tableCount(page, "chembl37.target_dictionary", 18552));
    await session.step(63, "When user opens the address \"/domains/chembl37/target_dictionary\"", () => openAddress(page, "/domains/chembl37/target_dictionary"), undefined, "Every ChEMBL target, read live, with the condition builder beside them");
    await session.step(64, "Then domain status count should contain text \"18,552\"", () => shouldContainText(page, el("domain status count"), "18,552"));
    await session.step(65, "When user enters \"Tyrosine-protein kinase BTK\" into domain search box", () => enterInto(page, "Tyrosine-protein kinase BTK", el("domain search box")));
    await session.step(66, "And user opens the \"Tyrosine-protein kinase BTK\" row of domain rows", () => openRow(page, "Tyrosine-protein kinase BTK"));
    await session.step(68, "Then domain breadcrumbs should contain text \"Tyrosine-protein kinase BTK\"", () => shouldContainText(page, el("domain breadcrumbs"), "Tyrosine-protein kinase BTK"), undefined, "BTK, and the assays measured on it, found through ChEMBL's own foreign keys");
    await session.step(69, "And domain children should contain text \"Assays (1662)\"", () => shouldContainText(page, el("domain children"), "Assays (1662)"));
    await session.step(71, "When user opens the address \"/domains/chembl37/molecule_dictionary\"", () => openAddress(page, "/domains/chembl37/molecule_dictionary"), undefined, "2.9 million compounds, searched live by name");
    await session.step(72, "And user enters \"IBRUTINIB\" into domain search box", () => enterInto(page, "IBRUTINIB", el("domain search box")));
    await session.step(73, "And user opens the \"CHEMBL1873475\" row of domain rows", () => openRow(page, "CHEMBL1873475"));
    await session.step(74, "Then domain breadcrumbs should contain text \"CHEMBL1873475\"", () => shouldContainText(page, el("domain breadcrumbs"), "CHEMBL1873475"));
    await session.step(76, "And the \"canonical_smiles\" column of the \"chembl37.compound_structures\" domain table should be recognized as \"Molecule\"", () => columnSemType(page, "canonical_smiles", "chembl37.compound_structures", "Molecule"), undefined, "Its structure lives in another table; Create Domain Schema recognized it as a molecule");
    await session.step(77, "When user clicks on \"Compound structures (1)\" tab in domain children", () => clickOn(page, el("\"Compound structures (1)\" tab in domain children")));
    await session.step(79, "And user clicks on \"Activities (1685)\" tab in domain children", () => clickOn(page, el("\"Activities (1685)\" tab in domain children")), undefined, "Its 1,685 activities, and one of them opened");
    await session.step(80, "And user clicks on the \"cell 1 of activity_id\" area of grid in domain children", () => clickArea(page, "cell 1 of activity_id", el("grid in domain children")));
    await session.step(81, "Then domain breadcrumbs should contain text \"6975312\"", () => shouldContainText(page, el("domain breadcrumbs"), "6975312"));
    await session.step(83, "When user follows the \"CHEMBL1794553\" link in domain form", () => followLink(page, "CHEMBL1794553"), undefined, "Every reference is a link: the assay it was measured in");
    await session.step(84, "Then domain breadcrumbs should contain text \"CHEMBL1794553\"", () => shouldContainText(page, el("domain breadcrumbs"), "CHEMBL1794553"));
    await session.step(85, "And domain form should contain text \"qHTS for Inhibitors of TGF-b\"", () => shouldContainText(page, el("domain form"), "qHTS for Inhibitors of TGF-b"));
    await session.step(87, "When user goes back", () => goBack(page), undefined, "…and back to the compound, one click away");
    await session.step(88, "And user follows the \"CHEMBL1201862\" link in domain form", () => followLink(page, "CHEMBL1201862"));
    await session.step(89, "Then domain breadcrumbs should contain text \"CHEMBL1201862\"", () => shouldContainText(page, el("domain breadcrumbs"), "CHEMBL1201862"));
    await session.step(90, "And domain form should contain text \"PubChem BioAssay data set\"", () => shouldContainText(page, el("domain form"), "PubChem BioAssay data set"));
    await session.step(92, "And domain children should contain text \"Activities (7434992)\"", () => shouldContainText(page, el("domain children"), "Activities (7434992)"), undefined, "The document behind it: 7.4 million activities, counted live in ChEMBL");
    await session.step(93, "And domain children should show their rows", () => childRows(page));
  });
});
