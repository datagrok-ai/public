/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sticky-meta/database-meta.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clearField, clickOn, enterInto, followingShouldBe, isExpanded, shouldBe, shouldContainText, shouldHaveValue, shouldNotBe, shouldNotHaveValue, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, dbMetaCleared, dbMetaOnServer, standHasReachableConnection} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Database meta of a database schema, of a table and of its column", () => {
  const session = feature(test, "features/sticky-meta/database-meta.feature", import.meta.url);
  test("Database meta of a database schema, of a table and of its column", {tag: ["@journey", "@serial", "@sticky-meta"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And the stand has a reachable \"PostgresTest\" connection", () => standHasReachableConnection(page, "PostgresTest"));
    await session.step(33, "And the Database meta of the \"PostgresTest\" connection is cleared now and at feature end", () => dbMetaCleared(page, "PostgresTest"));
    await session.step(34, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(36, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(37, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(38, "And Databases---Postgres---NorthwindTest tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest tree node inside browse tree")));
    await session.step(39, "And Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas tree node inside browse tree")));
    await run.scenario("A schema's Database meta is saved and cleared", async () => {
      await session.step(42, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
      await session.step(43, "Then context panel should contain text \"public\"", () => shouldContainText(page, el("context panel"), "public"));
      await session.step(44, "Given \"Database meta\" pane in context panel is expanded", () => isExpanded(page, el("\"Database meta\" pane in context panel")));
      await session.step(45, "Then Comment input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), ""));
      await session.step(46, "And \"LLM Comment\" input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel"), ""));
      await session.step(47, "When user enters \"test@#$! {time}\" into Comment input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("test@#$! {time}"), el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(48, "And user enters \"test@#$! llm {time}\" into \"LLM Comment\" input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("test@#$! llm {time}"), el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(49, "Then Save button in \"Database meta\" pane in context panel should be enabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "enabled"));
      await session.step(50, "When user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(51, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(52, "And the Database meta of \"public\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public", "PostgresTest", [["Comment",session.text("test@#$! {time}")],["LLM Comment",session.text("test@#$! llm {time}")]]), [["Comment",session.text("test@#$! {time}")],["LLM Comment",session.text("test@#$! llm {time}")]]);
      await session.step(55, "When user clears Comment input in \"Database meta\" pane in context panel", () => clearField(page, el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(56, "And user clears \"LLM Comment\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(57, "And user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(58, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(59, "And the Database meta of \"public\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public", "PostgresTest", [["Comment",""],["LLM Comment",""]]), [["Comment",""],["LLM Comment",""]]);
      await session.step(62, "And no errors should have been logged", () => noErrors(page));
      await session.step(63, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A table's Database meta is saved, and on no other table (4.1)", async () => {
      await session.step(66, "Given Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public tree node inside browse tree")));
      await session.step(67, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree")));
      await session.step(68, "Then the context panel should show \"categories\"", () => contextPanelShows(page, "categories"));
      await session.step(69, "Given \"Database meta\" pane in context panel is expanded", () => isExpanded(page, el("\"Database meta\" pane in context panel")));
      await session.step(70, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["Domains input in \"Database meta\" pane in context panel"],["\"Row Count\" input in \"Database meta\" pane in context panel"],["Comment input in \"Database meta\" pane in context panel"],["\"LLM Comment\" input in \"Database meta\" pane in context panel"],["Save button in \"Database meta\" pane in context panel"]]), [["Domains input in \"Database meta\" pane in context panel"],["\"Row Count\" input in \"Database meta\" pane in context panel"],["Comment input in \"Database meta\" pane in context panel"],["\"LLM Comment\" input in \"Database meta\" pane in context panel"],["Save button in \"Database meta\" pane in context panel"]]);
      await session.step(76, "And Comment input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), ""));
      await session.step(77, "And \"LLM Comment\" input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel"), ""));
      await session.step(78, "And \"Row Count\" input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("\"Row Count\" input in \"Database meta\" pane in context panel"), ""));
      await session.step(79, "When user enters \"bdd-sm-db-{time} table\" into Comment input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("bdd-sm-db-{time} table"), el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(80, "And user enters \"bdd-sm-db-{time} table llm\" into \"LLM Comment\" input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("bdd-sm-db-{time} table llm"), el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(81, "And user enters \"8\" into \"Row Count\" input in \"Database meta\" pane in context panel", () => enterInto(page, "8", el("\"Row Count\" input in \"Database meta\" pane in context panel")));
      await session.step(82, "Then Save button in \"Database meta\" pane in context panel should be enabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "enabled"));
      await session.step(83, "When user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(84, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(85, "And the Database meta of \"public.categories\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public.categories", "PostgresTest", [["Comment",session.text("bdd-sm-db-{time} table")],["LLM Comment",session.text("bdd-sm-db-{time} table llm")],["Row Count","8"]]), [["Comment",session.text("bdd-sm-db-{time} table")],["LLM Comment",session.text("bdd-sm-db-{time} table llm")],["Row Count","8"]]);
      await session.step(89, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---customers tree node inside browse tree")));
      await session.step(90, "Then the context panel should show \"customers\"", () => contextPanelShows(page, "customers"));
      await session.step(91, "And Comment input in \"Database meta\" pane in context panel should not have value \"bdd-sm-db-{time} table\"", () => shouldNotHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} table")));
      await session.step(92, "And \"LLM Comment\" input in \"Database meta\" pane in context panel should not have value \"bdd-sm-db-{time} table llm\"", () => shouldNotHaveValue(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} table llm")));
      await session.step(93, "And no errors should have been logged", () => noErrors(page));
      await session.step(94, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A column's Database meta is saved, and on no other column (4.2)", async () => {
      await session.step(97, "Given Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree")));
      await session.step(98, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree")));
      await session.step(99, "Then context panel should contain text \"categoryid\"", () => shouldContainText(page, el("context panel"), "categoryid"));
      await session.step(100, "Given \"Database meta\" pane in context panel is expanded", () => isExpanded(page, el("\"Database meta\" pane in context panel")));
      await session.step(101, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Is Unique\" input in \"Database meta\" pane in context panel"],["Min input in \"Database meta\" pane in context panel"],["Max input in \"Database meta\" pane in context panel"],["Values input in \"Database meta\" pane in context panel"],["\"Sample Values\" input in \"Database meta\" pane in context panel"],["\"Unique Count\" input in \"Database meta\" pane in context panel"],["Quality input in \"Database meta\" pane in context panel"],["Comment input in \"Database meta\" pane in context panel"],["\"LLM Comment\" input in \"Database meta\" pane in context panel"]]), [["\"Is Unique\" input in \"Database meta\" pane in context panel"],["Min input in \"Database meta\" pane in context panel"],["Max input in \"Database meta\" pane in context panel"],["Values input in \"Database meta\" pane in context panel"],["\"Sample Values\" input in \"Database meta\" pane in context panel"],["\"Unique Count\" input in \"Database meta\" pane in context panel"],["Quality input in \"Database meta\" pane in context panel"],["Comment input in \"Database meta\" pane in context panel"],["\"LLM Comment\" input in \"Database meta\" pane in context panel"]]);
      await session.step(111, "And \"Is Unique\" input in \"Database meta\" pane in context panel should not be checked", () => shouldNotBe(page, el("\"Is Unique\" input in \"Database meta\" pane in context panel"), "checked"));
      await session.step(112, "And Min input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("Min input in \"Database meta\" pane in context panel"), ""));
      await session.step(113, "And Comment input in \"Database meta\" pane in context panel should have value \"\"", () => shouldHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), ""));
      await session.step(114, "When user checks \"Is Unique\" input in \"Database meta\" pane in context panel", () => check(page, el("\"Is Unique\" input in \"Database meta\" pane in context panel")));
      await session.step(115, "And user enters \"1\" into Min input in \"Database meta\" pane in context panel", () => enterInto(page, "1", el("Min input in \"Database meta\" pane in context panel")));
      await session.step(116, "And user enters \"8\" into Max input in \"Database meta\" pane in context panel", () => enterInto(page, "8", el("Max input in \"Database meta\" pane in context panel")));
      await session.step(117, "And user enters \"1\" into Values input in \"Database meta\" pane in context panel", () => enterInto(page, "1", el("Values input in \"Database meta\" pane in context panel")));
      await session.step(118, "And user enters \"2\" into \"Sample Values\" input in \"Database meta\" pane in context panel", () => enterInto(page, "2", el("\"Sample Values\" input in \"Database meta\" pane in context panel")));
      await session.step(119, "And user enters \"8\" into \"Unique Count\" input in \"Database meta\" pane in context panel", () => enterInto(page, "8", el("\"Unique Count\" input in \"Database meta\" pane in context panel")));
      await session.step(120, "And user enters \"good {time}\" into Quality input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("good {time}"), el("Quality input in \"Database meta\" pane in context panel")));
      await session.step(121, "And user enters \"bdd-sm-db-{time} column\" into Comment input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("bdd-sm-db-{time} column"), el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(122, "And user enters \"bdd-sm-db-{time} column llm\" into \"LLM Comment\" input in \"Database meta\" pane in context panel", () => enterInto(page, session.text("bdd-sm-db-{time} column llm"), el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(123, "Then Save button in \"Database meta\" pane in context panel should be enabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "enabled"));
      await session.step(124, "When user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(125, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(126, "And the Database meta of \"public.categories.categoryid\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public.categories.categoryid", "PostgresTest", [["Is Unique","true"],["Min","1"],["Max","8"],["Values","1"],["Sample Values","2"],["Unique Count","8"],["Quality",session.text("good {time}")],["Comment",session.text("bdd-sm-db-{time} column")],["LLM Comment",session.text("bdd-sm-db-{time} column llm")]]), [["Is Unique","true"],["Min","1"],["Max","8"],["Values","1"],["Sample Values","2"],["Unique Count","8"],["Quality",session.text("good {time}")],["Comment",session.text("bdd-sm-db-{time} column")],["LLM Comment",session.text("bdd-sm-db-{time} column llm")]]);
      await session.step(136, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryname tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryname tree node inside browse tree")));
      await session.step(137, "Then context panel should contain text \"categoryname\"", () => shouldContainText(page, el("context panel"), "categoryname"));
      await session.step(138, "And Quality input in \"Database meta\" pane in context panel should not have value \"good {time}\"", () => shouldNotHaveValue(page, el("Quality input in \"Database meta\" pane in context panel"), session.text("good {time}")));
      await session.step(139, "And Comment input in \"Database meta\" pane in context panel should not have value \"bdd-sm-db-{time} column\"", () => shouldNotHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} column")));
      await session.step(140, "Given Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---NorthwindTest---Schemas---public---products tree node inside browse tree")));
      await session.step(141, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---products---categoryid tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---products---categoryid tree node inside browse tree")));
      await session.step(142, "Then context panel should contain text \"categoryid\"", () => shouldContainText(page, el("context panel"), "categoryid"));
      await session.step(143, "And Quality input in \"Database meta\" pane in context panel should not have value \"good {time}\"", () => shouldNotHaveValue(page, el("Quality input in \"Database meta\" pane in context panel"), session.text("good {time}")));
      await session.step(144, "And Comment input in \"Database meta\" pane in context panel should not have value \"bdd-sm-db-{time} column\"", () => shouldNotHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} column")));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The column's and the table's Database meta are shown again, cleared and saved (4.1, 4.2)", async () => {
      await session.step(148, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories---categoryid tree node inside browse tree")));
      await session.step(149, "Then context panel should contain text \"categoryid\"", () => shouldContainText(page, el("context panel"), "categoryid"));
      await session.step(150, "And \"Is Unique\" input in \"Database meta\" pane in context panel should be checked", () => shouldBe(page, el("\"Is Unique\" input in \"Database meta\" pane in context panel"), "checked"));
      await session.step(151, "And Comment input in \"Database meta\" pane in context panel should have value \"bdd-sm-db-{time} column\"", () => shouldHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} column")));
      await session.step(152, "When user unchecks \"Is Unique\" input in \"Database meta\" pane in context panel", () => uncheck(page, el("\"Is Unique\" input in \"Database meta\" pane in context panel")));
      await session.step(153, "And user clears Min input in \"Database meta\" pane in context panel", () => clearField(page, el("Min input in \"Database meta\" pane in context panel")));
      await session.step(154, "And user clears Max input in \"Database meta\" pane in context panel", () => clearField(page, el("Max input in \"Database meta\" pane in context panel")));
      await session.step(155, "And user clears Values input in \"Database meta\" pane in context panel", () => clearField(page, el("Values input in \"Database meta\" pane in context panel")));
      await session.step(156, "And user clears \"Sample Values\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"Sample Values\" input in \"Database meta\" pane in context panel")));
      await session.step(157, "And user clears \"Unique Count\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"Unique Count\" input in \"Database meta\" pane in context panel")));
      await session.step(158, "And user clears Quality input in \"Database meta\" pane in context panel", () => clearField(page, el("Quality input in \"Database meta\" pane in context panel")));
      await session.step(159, "And user clears Comment input in \"Database meta\" pane in context panel", () => clearField(page, el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(160, "And user clears \"LLM Comment\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(161, "And user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(162, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(163, "And the Database meta of \"public.categories.categoryid\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public.categories.categoryid", "PostgresTest", [["Min",""],["Max",""],["Values",""],["Unique Count",""],["Quality",""],["Comment",""],["LLM Comment",""]]), [["Min",""],["Max",""],["Values",""],["Unique Count",""],["Quality",""],["Comment",""],["LLM Comment",""]]);
      await session.step(171, "When user clicks on Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree", () => clickOn(page, el("Databases---Postgres---NorthwindTest---Schemas---public---categories tree node inside browse tree")));
      await session.step(172, "Then the context panel should show \"categories\"", () => contextPanelShows(page, "categories"));
      await session.step(173, "And Comment input in \"Database meta\" pane in context panel should have value \"bdd-sm-db-{time} table\"", () => shouldHaveValue(page, el("Comment input in \"Database meta\" pane in context panel"), session.text("bdd-sm-db-{time} table")));
      await session.step(174, "When user clears Comment input in \"Database meta\" pane in context panel", () => clearField(page, el("Comment input in \"Database meta\" pane in context panel")));
      await session.step(175, "And user clears \"LLM Comment\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"LLM Comment\" input in \"Database meta\" pane in context panel")));
      await session.step(176, "And user clears \"Row Count\" input in \"Database meta\" pane in context panel", () => clearField(page, el("\"Row Count\" input in \"Database meta\" pane in context panel")));
      await session.step(177, "And user clicks on Save button in \"Database meta\" pane in context panel", () => clickOn(page, el("Save button in \"Database meta\" pane in context panel")));
      await session.step(178, "Then Save button in \"Database meta\" pane in context panel should be disabled", () => shouldBe(page, el("Save button in \"Database meta\" pane in context panel"), "disabled"));
      await session.step(179, "And the Database meta of \"public.categories\" in the \"PostgresTest\" connection should be:", () => dbMetaOnServer(page, "public.categories", "PostgresTest", [["Comment",""],["LLM Comment",""],["Row Count",""]]), [["Comment",""],["LLM Comment",""],["Row Count",""]]);
      await session.step(183, "And no errors should have been logged", () => noErrors(page));
      await session.step(184, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
