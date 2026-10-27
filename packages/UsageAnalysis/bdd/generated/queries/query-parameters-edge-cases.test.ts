/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/queries/query-parameters-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.queries, core#229]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, enterInto, isExpanded, shouldBe, shouldNotOfferChoice, shouldOfferChoice} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, currentViewType, dialogCloses, queryTextOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Query parameters edge cases: the run dialog refuses what it cannot run", () => {
  const session = feature(test, "features/queries/query-parameters-edge-cases.feature", import.meta.url);
  test("An emptied integer and an emptied string both keep OK enabled", {tag: ["@edge-cases", "@realizes:views.queries", "@known-failure:42", "@realizes:core#229"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await knownFailure(async () => {
      await session.step(38, "When user clears Lim input in \"BDDQParamEdge{time}\" dialog", () => clearField(page, el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
      await session.step(39, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
      await session.step(40, "When user enters \"5\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "5", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
      await session.step(41, "And user clears \"Type Name\" input in \"BDDQParamEdge{time}\" dialog", () => clearField(page, el(session.text("\"Type Name\" input in \"BDDQParamEdge{time}\" dialog"))));
      await session.step(42, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
      await session.step(43, "When user clicks on CANCEL button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQParamEdge{time}\" dialog"))));
      await session.step(44, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
      await session.step(45, "And no errors should have been logged", () => noErrors(page));
    }, 42);
  });
  test("An integer above its maximum disables OK, and one inside the bounds enables it", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(48, "When user enters \"11\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "11", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(49, "Then OK button in \"BDDQParamEdge{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "disabled"));
    await session.step(50, "When user enters \"10\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "10", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(51, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(52, "When user enters \"0\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "0", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(53, "Then OK button in \"BDDQParamEdge{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "disabled"));
    await session.step(54, "When user clicks on CANCEL button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(55, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
  test("Text in an integer parameter disables OK", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(59, "When user enters \"abc\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "abc", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(60, "Then OK button in \"BDDQParamEdge{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "disabled"));
    await session.step(61, "When user enters \"3\" into Lim input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "3", el(session.text("Lim input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(62, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(63, "When user clicks on CANCEL button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(64, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
  });
  test("A value its validator rejects disables OK", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(68, "When user enters \"99\" into Code input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "99", el(session.text("Code input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(69, "Then OK button in \"BDDQParamEdge{time}\" dialog should be disabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "disabled"));
    await session.step(70, "When user enters \"1299\" into Code input in \"BDDQParamEdge{time}\" dialog", () => enterInto(page, "1299", el(session.text("Code input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(71, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(72, "When user clicks on CANCEL button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(73, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
    await session.step(74, "And no errors should have been logged", () => noErrors(page));
  });
  test("A choice parameter offers its choices only", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(77, "Then Kind input in \"BDDQParamEdge{time}\" dialog should offer the choice \"b\"", () => shouldOfferChoice(page, el(session.text("Kind input in \"BDDQParamEdge{time}\" dialog")), "b"));
    await session.step(78, "And Kind input in \"BDDQParamEdge{time}\" dialog should not offer the choice \"c\"", () => shouldNotOfferChoice(page, el(session.text("Kind input in \"BDDQParamEdge{time}\" dialog")), "c"));
    await session.step(79, "When user clicks on CANCEL button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(80, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
  });
  test("A nullable parameter left empty keeps OK enabled, and the query runs", {tag: ["@edge-cases", "@realizes:views.queries"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(18, "And a query \"BDDQParamEdge{time}\" on \"System:Datagrok\" is on the server:", () => queryTextOnServer(page, session.text("BDDQParamEdge{time}"), "System:Datagrok", "--input: string typeName = \"Project\"\n--input: int lim = 5 {min: 1; max: 10}\n--input: string kind = \"a\" {choices: [\"a\", \"b\"]}\n--input: string code = \"1234\" {validator: startsWith(value, \"12\")}\n--input: string note = \"x\" {nullable: true}\nselect name from entity_types where name = @typeName limit @lim"));
    await session.step(27, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(28, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(29, "And Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(30, "When user picks \"Run\" from the context menu of Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree", () => pickFromContextMenu(page, "Run", el(session.text("Databases---Postgres---Datagrok---BDDQParamEdge{time} tree node inside browse tree"))));
    await session.step(31, "Then \"BDDQParamEdge{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"BDDQParamEdge{time}\" dialog")), "visible"));
    await session.step(32, "And OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(84, "When user clears Note input in \"BDDQParamEdge{time}\" dialog", () => clearField(page, el(session.text("Note input in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(85, "Then OK button in \"BDDQParamEdge{time}\" dialog should be enabled", () => shouldBe(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog")), "enabled"));
    await session.step(86, "When user clicks on OK button in \"BDDQParamEdge{time}\" dialog", () => clickOn(page, el(session.text("OK button in \"BDDQParamEdge{time}\" dialog"))));
    await session.step(87, "Then the \"BDDQParamEdge{time}\" dialog should close", () => dialogCloses(page, session.text("BDDQParamEdge{time}")));
    await session.step(88, "And the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(89, "And the table should have 1 row", () => rowCount(page, 1));
    await session.step(90, "And no errors should have been logged", () => noErrors(page));
    await session.step(91, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
