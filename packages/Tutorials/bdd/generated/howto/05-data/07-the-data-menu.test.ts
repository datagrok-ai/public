/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/07-the-data-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-databases.js';
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
import {pointAtTopMenuItem} from '../../../bindings/howto-data.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/07-the-data-menu.feature", import.meta.url);
  test("What is in the Data menu", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "When the video says \"Data: commands that reshape and combine tables\" at the top", () => videoSaysAt(page, "Data: commands that reshape and combine tables", "top"));
    await session.step(21, "And user points at \"Data > Aggregate Rows...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Aggregate Rows..."), undefined, "Aggregate Rows: summarize groups of rows");
    await session.step(23, "And user points at \"Data > Unpivot...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Unpivot..."), undefined, "Unpivot: turn columns into rows");
    await session.step(25, "And user points at \"Data > Join Tables...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Join Tables..."), undefined, "Join, Link, Append, Compare: these need a second open table");
    await session.step(26, "Then \"Join Tables...\" menu item should be disabled", () => shouldBe(page, el("\"Join Tables...\" menu item"), "disabled"));
    await session.step(28, "When user points at \"Data > Batch Edit...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Batch Edit..."), undefined, "Batch Edit: one value into many rows at once");
    await session.step(30, "And user points at \"Data > Split Column...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Split Column..."), undefined, "Split Column: text into several columns, by a separator");
    await session.step(32, "And user points at \"Data > Split Column by RegExp...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Split Column by RegExp..."), undefined, "Split Column by RegExp: pull out parts that match a pattern");
    await session.step(34, "And user points at \"Data > Anonymize...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Anonymize..."), undefined, "Anonymize: replace values with made-up ones");
    await session.step(36, "And user points at \"Data > Categorize...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Categorize..."), undefined, "Categorize: bin numbers into named ranges");
    await session.step(38, "And user points at \"Data > Compare Columns...\" in the top menu", () => pointAtTopMenuItem(page, "Data > Compare Columns..."), undefined, "Compare Columns: row-by-row differences of two columns");
    await session.step(39, "Then \"Compare Columns...\" menu item should be enabled", () => shouldBe(page, el("\"Compare Columns...\" menu item"), "enabled"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
});
