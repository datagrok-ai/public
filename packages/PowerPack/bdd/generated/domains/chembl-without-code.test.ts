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
import {childRows, contextPanelClosed, followLink, headerPicture, openCard, openRow, pickAction, pickLink, pickTable, readDatabase, schemaGone, settingOn, synonymAbsent, synonymGone, synonymReads, tableCount} from '../../bindings/domains.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, enterInto, expand, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, dialogCloses, ownConnectionGone, shellSettingPutBack, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Connect ChEMBL and get an application over it without code", () => {
  const session = feature(test, "features/domains/chembl-without-code.feature", import.meta.url);
  test("Connect, name three tables, browse ChEMBL and change it from Datagrok", {tag: ["@guide", "@full-stand"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(14, "And the user's own connection \"ChEMBL37\" is removed now and at feature end", () => ownConnectionGone(page, "ChEMBL37"));
    await session.step(15, "And the domain schema \"chembl37\" is removed now and at feature end", () => schemaGone(page, "chembl37"));
    await session.step(16, "And the \"enableDomainDatabases\" shell setting is put back at feature end", () => shellSettingPutBack(page, "enableDomainDatabases"));
    await session.step(17, "And the \"enableDomainDatabases\" shell setting is on", () => settingOn(page, "enableDomainDatabases"));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(21, "When user picks \"New connection...\" from the context menu of Databases---Postgres tree node inside browse tree", () => pickFromContextMenu(page, "New connection...", el("Databases---Postgres tree node inside browse tree")), undefined, "The customer's ChEMBL 37: a Postgres database Datagrok has never seen");
    await session.step(22, "Then \"Add new connection\" dialog should be visible", () => shouldBe(page, el("\"Add new connection\" dialog"), "visible"));
    await session.step(23, "When user enters \"ChEMBL37\" into Name input in \"Add new connection\" dialog", () => enterInto(page, "ChEMBL37", el("Name input in \"Add new connection\" dialog")));
    await session.step(24, "And user enters \"host.docker.internal\" into Server input in \"Add new connection\" dialog", () => enterInto(page, "host.docker.internal", el("Server input in \"Add new connection\" dialog")));
    await session.step(25, "And user enters \"5441\" into Port input in \"Add new connection\" dialog", () => enterInto(page, "5441", el("Port input in \"Add new connection\" dialog")));
    await session.step(26, "And user enters \"chembl\" into Db input in \"Add new connection\" dialog", () => enterInto(page, "chembl", el("Db input in \"Add new connection\" dialog")));
    await session.step(27, "And user enters \"chembl_app\" into Login input in \"Add new connection\" dialog", () => enterInto(page, "chembl_app", el("Login input in \"Add new connection\" dialog")));
    await session.step(28, "And user enters \"chembl_app\" into Password input in \"Add new connection\" dialog", () => enterInto(page, "chembl_app", el("Password input in \"Add new connection\" dialog")));
    await session.step(29, "And user clicks on OK button in \"Add new connection\" dialog", () => clickOn(page, el("OK button in \"Add new connection\" dialog")));
    await session.step(30, "Then the \"Add new connection\" dialog should close", () => dialogCloses(page, "Add new connection"));
    await session.step(31, "And the synonym 9000001 is removed from the \"ChEMBL37\" database now and at feature end", () => synonymGone(page, 9000001, "ChEMBL37"));
    await session.step(33, "When user expands Databases---Postgres tree node inside browse tree", () => expand(page, el("Databases---Postgres tree node inside browse tree")), undefined, "Right-click its schema › Create Domain Schema…: every table, key and link is read from the database");
    await session.step(34, "And user expands Databases---Postgres---ChEMBL37 tree node inside browse tree", () => expand(page, el("Databases---Postgres---ChEMBL37 tree node inside browse tree")));
    await session.step(35, "And user expands Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree", () => expand(page, el("Databases---Postgres---ChEMBL37---Schemas tree node inside browse tree")));
    await session.step(36, "And user picks \"Create Domain Schema...\" from the context menu of Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree", () => pickFromContextMenu(page, "Create Domain Schema...", el("Databases---Postgres---ChEMBL37---Schemas---public tree node inside browse tree")));
    await session.step(37, "Then \"Create domain schema\" dialog should be visible", () => shouldBe(page, el("\"Create domain schema\" dialog"), "visible"));
    await session.step(38, "And \"Create domain schema\" dialog should contain text \"71 of 71 bindable tables\"", () => shouldContainText(page, el("\"Create domain schema\" dialog"), "71 of 71 bindable tables"));
    await session.step(40, "When user enters \"ChEMBL\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "ChEMBL", el("\"Friendly name\" input in \"Create domain schema\" dialog")), undefined, "Name it, let people write, and give three tables the names people use");
    await session.step(41, "And user enters \"chembl37\" into Identifier input in \"Create domain schema\" dialog", () => enterInto(page, "chembl37", el("Identifier input in \"Create domain schema\" dialog")));
    await session.step(42, "And user checks Writable input in \"Create domain schema\" dialog", () => check(page, el("Writable input in \"Create domain schema\" dialog")));
    await session.step(43, "And user clicks on the \"molecule_dictionary\" table in schema design", () => pickTable(page, "molecule_dictionary"));
    await session.step(44, "And user enters \"Compounds\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Compounds", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(45, "And user clicks on the \"target_dictionary\" table in schema design", () => pickTable(page, "target_dictionary"));
    await session.step(46, "And user enters \"Targets\" into \"Friendly name\" input in \"Create domain schema\" dialog", () => enterInto(page, "Targets", el("\"Friendly name\" input in \"Create domain schema\" dialog")));
    await session.step(47, "And user clicks on \"Next\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Next\" button in \"Create domain schema\" dialog")));
    await session.step(48, "And user clicks on \"Validate\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Validate\" button in \"Create domain schema\" dialog")));
    await session.step(49, "Then \"Create domain schema\" dialog should contain text \"No findings\"", () => shouldContainText(page, el("\"Create domain schema\" dialog"), "No findings"));
    await session.step(50, "When user clicks on \"Create\" button in \"Create domain schema\" dialog", () => clickOn(page, el("\"Create\" button in \"Create domain schema\" dialog")));
    await session.step(51, "Then the \"Create domain schema\" dialog should close", () => dialogCloses(page, "Create domain schema"));
    await session.step(52, "And the \"chembl37.molecule_dictionary\" domain table should hold 2921148 rows on the server", () => tableCount(page, "chembl37.molecule_dictionary", 2921148));
    await session.step(54, "And domain tables should contain text \"Compounds\"", () => shouldContainText(page, el("domain tables"), "Compounds"), undefined, "The schema's page: its 71 tables, registered under the names given them — nothing copied");
    await session.step(56, "When user double-clicks on the \"Compounds\" card of domain tables", () => openCard(page, "Compounds"), undefined, "2.9 million compounds, searched live in ChEMBL");
    await session.step(57, "Then domain status count should contain text \"2,921,148\"", () => shouldContainText(page, el("domain status count"), "2,921,148"));
    await session.step(58, "When user enters \"IBRUTINIB\" into domain search box", () => enterInto(page, "IBRUTINIB", el("domain search box")));
    await session.step(59, "And user opens the \"IBRUTINIB\" row of domain rows", () => openRow(page, "IBRUTINIB"));
    await session.step(61, "Then domain header should contain text \"IBRUTINIB\"", () => shouldContainText(page, el("domain header"), "IBRUTINIB"), undefined, "Its structure comes from compound_structures, found through ChEMBL's own keys");
    await session.step(62, "And domain header should show a structure", () => headerPicture(page));
    await session.step(63, "And domain children should contain text \"Activities (1685)\"", () => shouldContainText(page, el("domain children"), "Activities (1685)"));
    await session.step(64, "When user clicks on \"Molecule synonyms (20)\" tab in domain children", () => clickOn(page, el("\"Molecule synonyms (20)\" tab in domain children")));
    await session.step(65, "Then domain children should show their rows", () => childRows(page));
    await session.step(67, "When user clicks on \"Drug mechanism (1)\" tab in domain children", () => clickOn(page, el("\"Drug mechanism (1)\" tab in domain children")), undefined, "Every link is a page: its mechanism of action, and the target it acts on");
    await session.step(68, "Then domain children should show their rows", () => childRows(page));
    await session.step(69, "When user clicks on the \"cell 1 of mec_id\" area of grid in domain children", () => clickArea(page, "cell 1 of mec_id", el("grid in domain children")));
    await session.step(70, "Then domain breadcrumbs should contain text \"Drug mechanism\"", () => shouldContainText(page, el("domain breadcrumbs"), "Drug mechanism"));
    await session.step(71, "When user follows the \"Tyrosine-protein kinase BTK\" link in domain form", () => followLink(page, "Tyrosine-protein kinase BTK"));
    await session.step(72, "Then domain breadcrumbs should contain text \"Tyrosine-protein kinase BTK\"", () => shouldContainText(page, el("domain breadcrumbs"), "Tyrosine-protein kinase BTK"));
    await session.step(73, "And domain children should contain text \"Assays (1662)\"", () => shouldContainText(page, el("domain children"), "Assays (1662)"));
    await session.step(75, "When user opens the address \"/domains/chembl37/molecule_synonyms?entity=new\"", () => openAddress(page, "/domains/chembl37/molecule_synonyms?entity=new"), undefined, "A curator adds a trade name: a new row in ChEMBL's molecule_synonyms");
    await session.step(76, "And user picks \"IBRUTINIB\" in the \"Molregno\" link of domain form", () => pickLink(page, "IBRUTINIB", "Molregno"));
    await session.step(77, "And user enters \"TRADE_NAME\" into \"Syn type\" input in domain form", () => enterInto(page, "TRADE_NAME", el("\"Syn type\" input in domain form")));
    await session.step(78, "And user enters \"9000001\" into \"Molsyn id\" input in domain form", () => enterInto(page, "9000001", el("\"Molsyn id\" input in domain form")));
    await session.step(79, "And user enters \"IMBRUVICA (DEMO)\" into Synonyms input in domain form", () => enterInto(page, "IMBRUVICA (DEMO)", el("Synonyms input in domain form")));
    await session.step(80, "And user clicks on \"Save\" button in domain ribbon", () => clickOn(page, el("\"Save\" button in domain ribbon")));
    await session.step(81, "Then the synonym 9000001 in the \"ChEMBL37\" database should read \"IMBRUVICA (DEMO)\"", () => synonymReads(page, 9000001, "ChEMBL37", "IMBRUVICA (DEMO)"));
    await session.step(83, "When user reads \"select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001\" from the \"ChEMBL37\" database as \"In ChEMBL\"", () => readDatabase(page, "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001", "ChEMBL37", "In ChEMBL"), undefined, "…and it is in the customer's database");
    await session.step(84, "Then \"In ChEMBL\" dialog should be visible", () => shouldBe(page, el("\"In ChEMBL\" dialog"), "visible"));
    await session.step(85, "When user clicks on OK button in \"In ChEMBL\" dialog", () => clickOn(page, el("OK button in \"In ChEMBL\" dialog")));
    await session.step(87, "And user enters \"IMBRUVICA 560 MG (DEMO)\" into Synonyms input in domain form", () => enterInto(page, "IMBRUVICA 560 MG (DEMO)", el("Synonyms input in domain form")), undefined, "An edit, saved: ChEMBL holds the new text");
    await session.step(88, "And user clicks on \"Save\" button in domain ribbon", () => clickOn(page, el("\"Save\" button in domain ribbon")));
    await session.step(89, "Then the synonym 9000001 in the \"ChEMBL37\" database should read \"IMBRUVICA 560 MG (DEMO)\"", () => synonymReads(page, 9000001, "ChEMBL37", "IMBRUVICA 560 MG (DEMO)"));
    await session.step(90, "When user reads \"select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001\" from the \"ChEMBL37\" database as \"In ChEMBL after the edit\"", () => readDatabase(page, "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001", "ChEMBL37", "In ChEMBL after the edit"));
    await session.step(91, "Then \"In ChEMBL after the edit\" dialog should be visible", () => shouldBe(page, el("\"In ChEMBL after the edit\" dialog"), "visible"));
    await session.step(92, "When user clicks on OK button in \"In ChEMBL after the edit\" dialog", () => clickOn(page, el("OK button in \"In ChEMBL after the edit\" dialog")));
    await session.step(94, "And user picks \"Delete\" from domain actions", () => pickAction(page, "Delete"), undefined, "A delete, saved: the row is gone from ChEMBL");
    await session.step(95, "And user clicks on \"Save\" button in domain ribbon", () => clickOn(page, el("\"Save\" button in domain ribbon")));
    await session.step(96, "Then the synonym 9000001 should be gone from the \"ChEMBL37\" database", () => synonymAbsent(page, 9000001, "ChEMBL37"));
    await session.step(97, "When user reads \"select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001\" from the \"ChEMBL37\" database as \"In ChEMBL after the delete\"", () => readDatabase(page, "select molsyn_id, syn_type, synonyms from molecule_synonyms where molsyn_id = 9000001", "ChEMBL37", "In ChEMBL after the delete"));
    await session.step(98, "Then \"In ChEMBL after the delete\" dialog should be visible", () => shouldBe(page, el("\"In ChEMBL after the delete\" dialog"), "visible"));
    await session.step(99, "When user clicks on OK button in \"In ChEMBL after the delete\" dialog", () => clickOn(page, el("OK button in \"In ChEMBL after the delete\" dialog")));
  });
});
