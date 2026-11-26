/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/03-write-a-sql-query.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {runWithPlayIcon} from '../../../bindings/howto-databases.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, expand, hoverOver, replaceCode, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelClosed, currentViewType, noQueryOnServer, queriesOnServer, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/03-write-a-sql-query.feature", import.meta.url);
  test("Write a SQL query", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And no query named \"Assay types per category\" is on the server", () => noQueryOnServer(page, "Assay types per category"));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "When the video says \"Write your own SQL on a connection\"", () => videoSays(page, "Write your own SQL on a connection"));
    await session.step(22, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(23, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(25, "And user picks \"New Query...\" from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "New Query... opens the SQL editor for this connection");
    await session.step(26, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(28, "When user enters \"Assay types per category\" into Name input", () => enterInto(page, "Assay types per category", el("Name input")), undefined, "Give the query a name");
    await session.step(30, "And user replaces the code of code editor with \"select category, count(*) from biologics.assay_types group by category\"", () => replaceCode(page, el("code editor"), "select category, count(*) from biologics.assay_types group by category"), undefined, "Type the SQL");
    await session.step(32, "And user runs the query with the play icon", () => runWithPlayIcon(page), undefined, "Run it: the result shows under the editor");
    await session.step(33, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(34, "And the \"rows\" reading of grid should be 6", () => readingIs(page, "rows", el("grid"), 6));
    await session.step(36, "When user clicks on Save button", () => clickOn(page, el("Save button")), undefined, "Save keeps it on the server");
    await session.step(37, "Then 1 query named \"Assay types per category\" should be on the server", () => queriesOnServer(page, 1, "Assay types per category"));
    await session.step(38, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(40, "And user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "It is listed under its connection, ready to run");
    await session.step(41, "Then Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree should be visible", () => shouldBe(page, el("Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree"), "visible"));
    await session.step(42, "When user hovers over Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree", () => hoverOver(page, el("Databases---PostgresDart---Biologics---Assay-types-per-category tree node inside browse tree")));
    await session.step(43, "Then no errors should have been logged", () => noErrors(page));
  });
});
