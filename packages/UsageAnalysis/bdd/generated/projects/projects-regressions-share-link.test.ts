/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-share-link.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20930]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {onlyRowReads} from '../../bindings/projects-regressions.js';
import {datagrokQuery} from '../../bindings/projects-sources.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, noQueryOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: the share link of the Save project dialog", () => {
  const session = feature(test, "features/projects/projects-regressions-share-link.feature", import.meta.url);
  test("The share link follows the name typed into the Save dialog", {tag: ["@serial", "@realizes:views.projects", "@known-failure", "@realizes:GROK-20930"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(19, "And no query named \"BDDRegLinkQ{time}\" is on the server", () => noQueryOnServer(page, session.text("BDDRegLinkQ{time}")));
    await session.step(20, "And a query \"BDDRegLinkQ{time}\" on the Datagrok connection is:", () => datagrokQuery(page, session.text("BDDRegLinkQ{time}"), "--input: string typeName = \"Project\"\nselect name from entity_types where name = @typeName"));
    await session.step(25, "And Databases tree node inside browse tree is expanded", () => isExpanded(page, el("Databases tree node inside browse tree")));
    await session.step(26, "And Databases---Postgres tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres tree node inside browse tree")));
    await session.step(27, "When user clicks on \"Refresh\" icon inside browse toolbar", () => clickOn(page, el("\"Refresh\" icon inside browse toolbar")));
    await session.step(28, "Given Databases---Postgres---Datagrok tree node inside browse tree is expanded", () => isExpanded(page, el("Databases---Postgres---Datagrok tree node inside browse tree")));
    await session.step(29, "When user double-clicks Databases---Postgres---Datagrok---BDDRegLinkQ{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("Databases---Postgres---Datagrok---BDDRegLinkQ{time} tree node inside browse tree"))));
    await session.step(30, "Then the \"BDDRegLinkQ{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDRegLinkQ{time}")));
    await session.step(31, "And the only row of the table should read \"Project\" in the \"name\" column", () => onlyRowReads(page, "Project", "name"));
    await session.step(32, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(33, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(34, "And save dialog share link should contain text \"BDDRegLinkQ{time}?typeName=Project\"", () => shouldContainText(page, el("save dialog share link"), session.text("BDDRegLinkQ{time}?typeName=Project")));
    await session.step(35, "When user enters \"BDDRegShareName{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegShareName{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(36, "And user presses Tab", () => pressKey(page, "Tab"));
    await session.step(37, "Then Name text input in \"Save project\" dialog should have value \"BDDRegShareName{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("BDDRegShareName{time}")));
    await knownFailure(async () => {
      await session.step(41, "Then save dialog share link should contain text \"BDDRegShareName{time}?typeName=Project\"", () => shouldContainText(page, el("save dialog share link"), session.text("BDDRegShareName{time}?typeName=Project")));
    });
  });
});
