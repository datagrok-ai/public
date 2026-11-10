/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/identifiers/identifiers-catalogs.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.db-explorer]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {findCatalog, identifiersConfigCovers, noIdentifiersConfig, savedIdentifiers} from '../../bindings/identifiers.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue, shouldNotContainText, shouldNotOfferChoice, shouldOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, contextPanelOpen, openTableOf, standHasReachableConnection, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Identifiers of a database with catalogs", () => {
  const session = feature(test, "features/identifiers/identifiers-catalogs.feature", import.meta.url);
  test.skip("Identifiers of a database with catalogs", {tag: ["@skip", "@journey", "@serial", "@realizes:powerpack.db-explorer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(30, "Given user is logged in", () => loggedIn(page));
    await session.step(31, "And the stand has a reachable \"CatalogsDemo\" connection", () => standHasReachableConnection(page, "CatalogsDemo"));
    await session.step(32, "And no identifiers configuration of the \"CatalogsDemo\" connection is on the server", () => noIdentifiersConfig(page, "CatalogsDemo"));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The dialog asks for the catalog first and opens on the connection's own", async () => {
      await session.step(36, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(37, "And Databases---MS-SQL tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL tree node inside browse tree")));
      await session.step(38, "When user picks \"Configure Identifiers...\" from the context menu of Databases---MS-SQL---Catalogs-Demo tree node inside browse tree", () => pickFromContextMenu(page, "Configure Identifiers...", el("Databases---MS-SQL---Catalogs-Demo tree node inside browse tree")));
      await session.step(39, "Then \"Select catalog and schema for Identifiers Configuration\" dialog should be visible", () => shouldBe(page, el("\"Select catalog and schema for Identifiers Configuration\" dialog"), "visible"));
      await session.step(40, "And Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog should have value \"archive\"", () => shouldHaveValue(page, el("Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "archive"));
      await session.step(41, "And Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog should offer the choice \"chem\"", () => shouldOfferChoice(page, el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "chem"));
      await session.step(42, "And Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog should not offer the choice \"bio\"", () => shouldNotOfferChoice(page, el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "bio"));
      await session.step(43, "And \"Select catalog and schema for Identifiers Configuration\" dialog should contain text \"Select the schema that holds the identifier tables.\"", () => shouldContainText(page, el("\"Select catalog and schema for Identifiers Configuration\" dialog"), "Select the schema that holds the identifier tables."));
      await session.step(44, "And OK button in \"Select catalog and schema for Identifiers Configuration\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Select catalog and schema for Identifiers Configuration\" dialog"), "disabled"));
    });
    await run.scenario("Among two hundred catalogs one is found by a part of its name", async () => {
      await session.step(47, "Then \"Select catalog and schema for Identifiers Configuration\" dialog should contain text \"207 catalogs; archive is the connection's own.\"", () => shouldContainText(page, el("\"Select catalog and schema for Identifiers Configuration\" dialog"), "207 catalogs; archive is the connection's own."));
      await session.step(48, "When user clicks on catalog search icon", () => clickOn(page, el("catalog search icon")));
      await session.step(49, "And user finds the catalog \"tenant_150\" by typing \"tenant_15\"", () => findCatalog(page, "tenant_150", "tenant_15"));
      await session.step(50, "Then Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog should have value \"tenant_150\"", () => shouldHaveValue(page, el("Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "tenant_150"));
      await session.step(51, "And Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog should have value \"dbo\"", () => shouldHaveValue(page, el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "dbo"));
      await session.step(52, "And OK button in \"Select catalog and schema for Identifiers Configuration\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Select catalog and schema for Identifiers Configuration\" dialog"), "enabled"));
    });
    await run.scenario("The editor of the connection's own catalog names it and says nothing of links", async () => {
      await session.step(55, "When user selects \"archive\" in Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog", () => selectIn(page, "archive", el("Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(56, "Then Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog should offer the choice \"core\"", () => shouldOfferChoice(page, el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "core"));
      await session.step(57, "When user selects \"chem\" in Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog", () => selectIn(page, "chem", el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(58, "And user clicks on OK button in \"Select catalog and schema for Identifiers Configuration\" dialog", () => clickOn(page, el("OK button in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(59, "Then the \"CatalogsDemo.archive.chem Identifiers\" view should be current", () => viewIsCurrent(page, "CatalogsDemo.archive.chem Identifiers"));
      await session.step(60, "And identifiers editor should contain text \"Catalogs Demo: archive.chem\"", () => shouldContainText(page, el("identifiers editor"), "Catalogs Demo: archive.chem"));
      await session.step(61, "And identifiers editor should not contain text \"add them under Explicit References\"", () => shouldNotContainText(page, el("identifiers editor"), "add them under Explicit References"));
    });
    await run.scenario("Another catalog brings its own schemas, and its editor says where links come from", async () => {
      await session.step(64, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(65, "When user picks \"Configure Identifiers...\" from the context menu of Databases---MS-SQL---Catalogs-Demo tree node inside browse tree", () => pickFromContextMenu(page, "Configure Identifiers...", el("Databases---MS-SQL---Catalogs-Demo tree node inside browse tree")));
      await session.step(66, "And user selects \"research\" in Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog", () => selectIn(page, "research", el("Catalog input in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(67, "Then Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog should offer the choice \"bio\"", () => shouldOfferChoice(page, el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog"), "bio"));
      await session.step(68, "When user selects \"chem\" in Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog", () => selectIn(page, "chem", el("Schema input in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(69, "And user clicks on OK button in \"Select catalog and schema for Identifiers Configuration\" dialog", () => clickOn(page, el("OK button in \"Select catalog and schema for Identifiers Configuration\" dialog")));
      await session.step(70, "Then the \"CatalogsDemo.research.chem Identifiers\" view should be current", () => viewIsCurrent(page, "CatalogsDemo.research.chem Identifiers"));
      await session.step(71, "And identifiers editor should contain text \"Catalogs Demo: research.chem\"", () => shouldContainText(page, el("identifiers editor"), "Catalogs Demo: research.chem"));
      await session.step(72, "And identifiers editor should contain text \"For research, add them under Explicit References.\"", () => shouldContainText(page, el("identifiers editor"), "For research, add them under Explicit References."));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
      await session.step(74, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A configuration made for another catalog brings its identifiers and leaves the editor where it is", async () => {
      await session.step(77, "When user clicks on Import JSON button in identifiers editor", () => clickOn(page, el("Import JSON button in identifiers editor")));
      await session.step(78, "Then \"Import Configuration\" dialog should be visible", () => shouldBe(page, el("\"Import Configuration\" dialog"), "visible"));
      await session.step(79, "When user enters '{\"connectionName\":\"CatalogsDemo\",\"nqName\":\"CatalogsDemo\",\"dataSourceName\":\"MS SQL\",\"catalogName\":\"sales\",\"schemaName\":\"dbo\",\"entryPoints\":{\"BDD_IMPORTED_ID\":{\"table\":\"batches\",\"column\":\"batch_id\"}}}' into configuration text", () => enterInto(page, "{\"connectionName\":\"CatalogsDemo\",\"nqName\":\"CatalogsDemo\",\"dataSourceName\":\"MS SQL\",\"catalogName\":\"sales\",\"schemaName\":\"dbo\",\"entryPoints\":{\"BDD_IMPORTED_ID\":{\"table\":\"batches\",\"column\":\"batch_id\"}}}", el("configuration text")));
      await session.step(80, "And user clicks on Import button in \"Import Configuration\" dialog", () => clickOn(page, el("Import button in \"Import Configuration\" dialog")));
      await session.step(81, "Then an info balloon containing \"imported into CatalogsDemo.research.chem\" should have been shown", () => infoBalloonText(page, "imported into CatalogsDemo.research.chem"));
      await session.step(82, "And identifiers editor should contain text \"BDD_IMPORTED_ID\"", () => shouldContainText(page, el("identifiers editor"), "BDD_IMPORTED_ID"));
      await session.step(83, "And identifiers editor should contain text \"Catalogs Demo: research.chem\"", () => shouldContainText(page, el("identifiers editor"), "Catalogs Demo: research.chem"));
    });
    await run.scenario("An identifier is previewed from the picked catalog and saved with it", async () => {
      await session.step(86, "When user clicks on add identifier icon", () => clickOn(page, el("add identifier icon")));
      await session.step(87, "Then \"Add Identifier\" dialog should be visible", () => shouldBe(page, el("\"Add Identifier\" dialog"), "visible"));
      await session.step(88, "And Table input in \"Add Identifier\" dialog should offer the choice \"batches\"", () => shouldOfferChoice(page, el("Table input in \"Add Identifier\" dialog"), "batches"));
      await session.step(89, "When user enters \"BDD_COMPOUND_ID\" into Semantic Type input in \"Add Identifier\" dialog", () => enterInto(page, "BDD_COMPOUND_ID", el("Semantic Type input in \"Add Identifier\" dialog")));
      await session.step(90, "And user selects \"compounds\" in Table input in \"Add Identifier\" dialog", () => selectIn(page, "compounds", el("Table input in \"Add Identifier\" dialog")));
      await session.step(91, "And user selects \"compound_id\" in Column input in \"Add Identifier\" dialog", () => selectIn(page, "compound_id", el("Column input in \"Add Identifier\" dialog")));
      await session.step(92, "Then identifier preview should contain text \"DGC-\"", () => shouldContainText(page, el("identifier preview"), "DGC-"));
      await session.step(93, "And identifier preview should not contain text \"ARCHIVED\"", () => shouldNotContainText(page, el("identifier preview"), "ARCHIVED"));
      await session.step(94, "When user enters \"DGC-\\d{6}\" into Match Regexp input in \"Add Identifier\" dialog", () => enterInto(page, "DGC-\\d{6}", el("Match Regexp input in \"Add Identifier\" dialog")));
      await session.step(95, "And user clicks on ADD button in \"Add Identifier\" dialog", () => clickOn(page, el("ADD button in \"Add Identifier\" dialog")));
      await session.step(96, "Then identifiers editor should contain text \"BDD_COMPOUND_ID\"", () => shouldContainText(page, el("identifiers editor"), "BDD_COMPOUND_ID"));
      await session.step(97, "When user clicks on Save button in identifiers editor", () => clickOn(page, el("Save button in identifiers editor")));
      await session.step(98, "Then an info balloon containing \"Identifier configuration saved\" should have been shown", () => infoBalloonText(page, "Identifier configuration saved"));
      await session.step(99, "And the identifiers configurations of the \"CatalogsDemo\" connection on the server should cover \"research.chem\"", () => identifiersConfigCovers(page, "CatalogsDemo", "research.chem"));
      await session.step(100, "And the saved identifiers of \"research.chem\" in the \"CatalogsDemo\" connection should be \"BDD_COMPOUND_ID, BDD_IMPORTED_ID\"", () => savedIdentifiers(page, "research.chem", "CatalogsDemo", "BDD_COMPOUND_ID, BDD_IMPORTED_ID"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The configuration is listed under the connection by catalog and schema", async () => {
      await session.step(104, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(105, "And Databases---MS-SQL---Catalogs-Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---MS-SQL---Catalogs-Demo tree node inside browse tree")));
      await session.step(106, "Then Databases---MS-SQL---Catalogs-Demo---Identifiers---research.chem tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---MS-SQL---Catalogs-Demo---Identifiers---research.chem tree node inside browse tree"), "visible"));
    });
    await run.scenario("After a reload the identifier opens the row of the picked catalog", async () => {
      await session.step(109, "When user reloads the page", () => reloadPage(page));
      await session.step(110, "Given the package autostarts have completed", () => autostartsCompleted(page));
      await session.step(111, "And the context panel is open", () => contextPanelOpen(page));
      await session.step(112, "And user opens a table \"bdd-compounds\" with:", () => openTableOf(page, "bdd-compounds", [["compound","note"],["DGC-000001","first"],["DGC-000005","second"]]), [["compound","note"],["DGC-000001","first"],["DGC-000005","second"]]);
      await session.step(116, "Then \"compound\" column should have semantic type \"BDD_COMPOUND_ID\"", () => columnSemType(page, "compound", "BDD_COMPOUND_ID"));
      await session.step(117, "When user clicks on the \"cell 1 of compound\" area of grid", () => clickArea(page, "cell 1 of compound", el("grid")));
      await session.step(118, "Then context panel should contain text \"Aspirin\"", () => shouldContainText(page, el("context panel"), "Aspirin"));
      await session.step(119, "And context panel should not contain text \"ARCHIVED\"", () => shouldNotContainText(page, el("context panel"), "ARCHIVED"));
      await session.step(120, "And no errors should have been logged", () => noErrors(page));
      await session.step(121, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
