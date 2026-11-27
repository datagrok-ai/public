/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/identifiers/identifiers.feature
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
import {identifiersConfigCovers, noIdentifiersConfig, savedEntityNames, savedIdentifiers} from '../../bindings/identifiers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Identifiers of a database connection", () => {
  const session = feature(test, "features/identifiers/identifiers.feature", import.meta.url);
  test.skip("Identifiers of a database connection", {tag: ["@skip", "@journey", "@serial", "@realizes:powerpack.db-explorer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(23, "Given user is logged in", () => loggedIn(page));
    await session.step(24, "And no identifiers configuration of the \"System:Datagrok\" connection is on the server", () => noIdentifiersConfig(page, "System:Datagrok"));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("The dialog asks for the schema alone and opens on public", async () => {
      await session.step(28, "Given Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
      await session.step(29, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
      await session.step(30, "When user picks \"Configure Identifiers...\" from the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => pickFromContextMenu(page, "Configure Identifiers...", el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(31, "Then \"Select primary schema for Identifiers Configuration\" dialog should be visible", () => shouldBe(page, el("\"Select primary schema for Identifiers Configuration\" dialog"), "visible"));
      await session.step(32, "And Schema input in \"Select primary schema for Identifiers Configuration\" dialog should have value \"public\"", () => shouldHaveValue(page, el("Schema input in \"Select primary schema for Identifiers Configuration\" dialog"), "public"));
      await session.step(33, "And Catalog input in \"Select primary schema for Identifiers Configuration\" dialog should be absent", () => shouldBe(page, el("Catalog input in \"Select primary schema for Identifiers Configuration\" dialog"), "absent"));
      await session.step(34, "When user clicks on OK button in \"Select primary schema for Identifiers Configuration\" dialog", () => clickOn(page, el("OK button in \"Select primary schema for Identifiers Configuration\" dialog")));
      await session.step(35, "Then the \"Datagrok.public Identifiers\" view should be current", () => viewIsCurrent(page, "Datagrok.public Identifiers"));
      await session.step(36, "And identifiers editor should be visible", () => shouldBe(page, el("identifiers editor"), "visible"));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
      await session.step(38, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("An identifier is previewed from its table and saved", async () => {
      await session.step(41, "When user clicks on add identifier icon", () => clickOn(page, el("add identifier icon")));
      await session.step(42, "Then \"Add Identifier\" dialog should be visible", () => shouldBe(page, el("\"Add Identifier\" dialog"), "visible"));
      await session.step(43, "When user enters \"BDD_EVENT_TYPE\" into Semantic Type input in \"Add Identifier\" dialog", () => enterInto(page, "BDD_EVENT_TYPE", el("Semantic Type input in \"Add Identifier\" dialog")));
      await session.step(44, "And user selects \"event_types\" in Table input in \"Add Identifier\" dialog", () => selectIn(page, "event_types", el("Table input in \"Add Identifier\" dialog")));
      await session.step(45, "And user selects \"friendly_name\" in Column input in \"Add Identifier\" dialog", () => selectIn(page, "friendly_name", el("Column input in \"Add Identifier\" dialog")));
      await session.step(46, "Then identifier preview should contain text \"friendly_name\"", () => shouldContainText(page, el("identifier preview"), "friendly_name"));
      await session.step(47, "When user clicks on ADD button in \"Add Identifier\" dialog", () => clickOn(page, el("ADD button in \"Add Identifier\" dialog")));
      await session.step(48, "Then identifiers editor should contain text \"BDD_EVENT_TYPE\"", () => shouldContainText(page, el("identifiers editor"), "BDD_EVENT_TYPE"));
      await session.step(49, "When user clicks on Save button in identifiers editor", () => clickOn(page, el("Save button in identifiers editor")));
      await session.step(50, "Then an info balloon containing \"Identifier configuration saved\" should have been shown", () => infoBalloonText(page, "Identifier configuration saved"));
      await session.step(51, "And the identifiers configurations of the \"System:Datagrok\" connection on the server should cover \"public\"", () => identifiersConfigCovers(page, "System:Datagrok", "public"));
      await session.step(52, "And the saved identifiers of \"public\" in the \"System:Datagrok\" connection should be \"BDD_EVENT_TYPE\"", () => savedIdentifiers(page, "public", "System:Datagrok", "BDD_EVENT_TYPE"));
      await session.step(53, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table is named for the Links section of the cards and saved with the identifiers", async () => {
      await session.step(58, "When user expands Entity Names section inside identifiers editor", () => expand(page, el("Entity Names section inside identifiers editor")));
      await session.step(59, "And user clicks on add entity name icon", () => clickOn(page, el("add entity name icon")));
      await session.step(60, "Then \"Add Entity Name\" dialog should be visible", () => shouldBe(page, el("\"Add Entity Name\" dialog"), "visible"));
      await session.step(61, "When user selects \"event_types\" in Table input in \"Add Entity Name\" dialog", () => selectIn(page, "event_types", el("Table input in \"Add Entity Name\" dialog")));
      await session.step(62, "And user enters \"Event type\" into Name input in \"Add Entity Name\" dialog", () => enterInto(page, "Event type", el("Name input in \"Add Entity Name\" dialog")));
      await session.step(63, "And user enters \"EVT\" into Code input in \"Add Entity Name\" dialog", () => enterInto(page, "EVT", el("Code input in \"Add Entity Name\" dialog")));
      await session.step(64, "And user clicks on ADD button in \"Add Entity Name\" dialog", () => clickOn(page, el("ADD button in \"Add Entity Name\" dialog")));
      await session.step(65, "Then identifiers editor should contain text \"Event type\"", () => shouldContainText(page, el("identifiers editor"), "Event type"));
      await session.step(67, "When user clicks on Save button in identifiers editor", () => clickOn(page, el("Save button in identifiers editor")));
      await session.step(68, "And user clicks on OK button in \"Overwrite existing configuration?\" dialog", () => clickOn(page, el("OK button in \"Overwrite existing configuration?\" dialog")));
      await session.step(69, "Then an info balloon containing \"Identifier configuration saved\" should have been shown", () => infoBalloonText(page, "Identifier configuration saved"));
      await session.step(70, "And the saved entity names of \"public\" in the \"System:Datagrok\" connection should be \"event_types: Event type (EVT)\"", () => savedEntityNames(page, "public", "System:Datagrok", "event_types: Event type (EVT)"));
      await session.step(71, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The saved configuration is listed under the connection and opens from there", async () => {
      await session.step(74, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(75, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(76, "When user clicks on Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree")));
      await session.step(77, "Then identifiers editor should contain text \"BDD_EVENT_TYPE\"", () => shouldContainText(page, el("identifiers editor"), "BDD_EVENT_TYPE"));
      await session.step(78, "When user expands Entity Names section inside identifiers editor", () => expand(page, el("Entity Names section inside identifiers editor")));
      await session.step(79, "Then identifiers editor should contain text \"Event type\"", () => shouldContainText(page, el("identifiers editor"), "Event type"));
      await session.step(80, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Configure Identifiers on the same schema opens the saved configuration", async () => {
      await session.step(83, "When user picks \"Configure Identifiers...\" from the context menu of Databases---Postgres---Datagrok tree node inside browse tree", () => pickFromContextMenu(page, "Configure Identifiers...", el("Databases---Postgres---Datagrok tree node inside browse tree")));
      await session.step(84, "Then \"Select primary schema for Identifiers Configuration\" dialog should contain text \"public already has identifiers. OK opens them for editing.\"", () => shouldContainText(page, el("\"Select primary schema for Identifiers Configuration\" dialog"), "public already has identifiers. OK opens them for editing."));
      await session.step(85, "When user clicks on OK button in \"Select primary schema for Identifiers Configuration\" dialog", () => clickOn(page, el("OK button in \"Select primary schema for Identifiers Configuration\" dialog")));
      await session.step(86, "Then an info balloon containing \"already exists\" should have been shown", () => infoBalloonText(page, "already exists"));
      await session.step(87, "And the \"Datagrok.public Identifiers\" view should be current", () => viewIsCurrent(page, "Datagrok.public Identifiers"));
      await session.step(88, "And identifiers editor should contain text \"BDD_EVENT_TYPE\"", () => shouldContainText(page, el("identifiers editor"), "BDD_EVENT_TYPE"));
    });
    await run.scenario("Remove takes the configuration off the connection and the server", async () => {
      await session.step(91, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(92, "When user picks \"Remove...\" from the context menu of Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree", () => pickFromContextMenu(page, "Remove...", el("Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree")));
      await session.step(93, "Then \"Remove Identifier Configuration\" dialog should be visible", () => shouldBe(page, el("\"Remove Identifier Configuration\" dialog"), "visible"));
      await session.step(94, "When user clicks on OK button in \"Remove Identifier Configuration\" dialog", () => clickOn(page, el("OK button in \"Remove Identifier Configuration\" dialog")));
      await session.step(95, "Then an info balloon containing \"Identifier configuration removed\" should have been shown", () => infoBalloonText(page, "Identifier configuration removed"));
      await session.step(96, "And Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree should be absent", () => shouldBe(page, el("Databases---Postgres---Datagrok---Identifiers---public tree node inside browse tree"), "absent"));
      await session.step(97, "And the identifiers configurations of the \"System:Datagrok\" connection on the server should cover \"\"", () => identifiersConfigCovers(page, "System:Datagrok", ""));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
