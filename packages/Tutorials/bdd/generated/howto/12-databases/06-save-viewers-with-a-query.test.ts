/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/12-databases/06-save-viewers-with-a-query.feature
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
import {ownLayoutsGone} from '../../../bindings/howto-dashboards.js';
import {runWithPlayIcon} from '../../../bindings/howto-databases.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, enterInto, expand, hoverOver, replaceCode, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, contextPanelClosed, currentViewType, entityHasLayout, noQueryOnServer, queriesOnServer, simpleModeOff, toolboxPaneHidden, tooltipHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Databases", () => {
  const session = feature(test, "features/howto/12-databases/06-save-viewers-with-a-query.feature", import.meta.url);
  test("Save viewers with a query", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And no query named \"Expression yields\" is on the server", () => noQueryOnServer(page, "Expression yields"));
    await session.step(19, "And the user's own layouts named \"Expression yields, ExpressionYields\" are removed now and at feature end", () => ownLayoutsGone(page, "Expression yields, ExpressionYields"));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "When the video says \"A query can keep its viewers: its result opens with them\"", () => videoSays(page, "A query can keep its viewers: its result opens with them"));
    await session.step(24, "And user expands Databases tree node inside browse tree", () => expand(page, el("Databases tree node inside browse tree")));
    await session.step(25, "And user expands Databases---PostgresDart tree node inside browse tree", () => expand(page, el("Databases---PostgresDart tree node inside browse tree")));
    await session.step(26, "And user picks \"New Query...\" from the context menu of Databases---PostgresDart---Biologics tree node inside browse tree", () => pickFromContextMenu(page, "New Query...", el("Databases---PostgresDart---Biologics tree node inside browse tree")));
    await session.step(27, "Then the current view should be a DataQueryView view", () => currentViewType(page, "DataQueryView"));
    await session.step(28, "When user enters \"Expression yields\" into Name input", () => enterInto(page, "Expression yields", el("Name input")));
    await session.step(29, "And user replaces the code of code editor with \"select expression_system, chain, yield_mg from biologics.expression_batches\"", () => replaceCode(page, el("code editor"), "select expression_system, chain, yield_mg from biologics.expression_batches"));
    await session.step(31, "And user clicks on Layout tab", () => clickOn(page, el("Layout tab")), undefined, "The Layout tab sets how the result looks");
    await session.step(32, "Then \"Run query to get data and edit layout\" text should be visible", () => shouldBe(page, el("\"Run query to get data and edit layout\" text"), "visible"));
    await session.step(34, "When user runs the query with the play icon", () => runWithPlayIcon(page), undefined, "Run the query to fill it");
    await session.step(35, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    await session.step(37, "When user clicks on box plot icon in toolbox", () => clickOn(page, el("box plot icon in toolbox")), undefined, "Add viewers from the toolbox");
    await session.step(38, "And user clicks on histogram icon in toolbox", () => clickOn(page, el("histogram icon in toolbox")));
    await session.step(39, "Then box plot viewer should be visible", () => shouldBe(page, el("box plot viewer"), "visible"));
    await session.step(40, "And histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
    await session.step(42, "When user clicks on Save button", () => clickOn(page, el("Save button")), undefined, "SAVE keeps them with the query");
    await session.step(43, "Then 1 query named \"Expression yields\" should be on the server", () => queriesOnServer(page, 1, "Expression yields"));
    await session.step(44, "And the query \"Expression yields\" on the server should have a layout", () => entityHasLayout(page, "query", "Expression yields"));
    await session.step(46, "When user closes \"Expression yields\" view tab", () => close(page, el("\"Expression yields\" view tab")), undefined, "Close the editor");
    await session.step(47, "Then \"Expression yields\" view tab should be absent", () => shouldBe(page, el("\"Expression yields\" view tab"), "absent"));
    await session.step(50, "Given the tooltip is hidden", () => tooltipHidden(page));
    await session.step(52, "When user expands Databases---PostgresDart---Biologics tree node inside browse tree", () => expand(page, el("Databases---PostgresDart---Biologics tree node inside browse tree")), undefined, "Run the saved query from Browse");
    await session.step(53, "And user picks \"Run\" from the context menu of Databases---PostgresDart---Biologics---Expression-yields tree node inside browse tree", () => pickFromContextMenu(page, "Run", el("Databases---PostgresDart---Biologics---Expression-yields tree node inside browse tree")));
    await session.step(54, "Then the current view should be a TableView view", () => currentViewType(page, "TableView"));
    await session.step(55, "And the table should have 1200 rows", () => rowCount(page, 1200));
    await session.step(56, "And the open tableview should have 1 box plot viewer", () => viewerCount(page, 1, "box plot"));
    await session.step(57, "And the open tableview should have 1 histogram viewer", () => viewerCount(page, 1, "histogram"));
    await session.step(58, "When the video says \"The result opens with the viewers saved with the query\" at the top", () => videoSaysAt(page, "The result opens with the viewers saved with the query", "top"));
    await session.step(59, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
    await session.step(60, "Then no errors should have been logged", () => noErrors(page));
  });
});
