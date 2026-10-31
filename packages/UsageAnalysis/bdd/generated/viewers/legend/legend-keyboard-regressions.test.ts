/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/legend/legend-keyboard-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.legend, public#2548, core#1159]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, legendLists, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Legend regressions — the arrow keys", () => {
  const session = feature(test, "features/viewers/legend/legend-keyboard-regressions.feature", import.meta.url);
  test("Down and Up move the exclusive selection to the next and the previous category (public#2548)", {tag: ["@viewers", "@realizes:viewers.legend", "@known-failure:26", "@realizes:public#2548", "@realizes:core#1159"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]), [["X","WEIGHT"],["Y","HEIGHT"],["Color","RACE"],["Legend Visibility","Always"],["Legend Position","Right"]]);
    await session.step(18, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
    await knownFailure(async () => {
      await session.step(22, "When user clicks on \"Asian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Asian", el("scatter plot viewer")));
      await session.step(23, "Then \"Asian\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(24, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(25, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
      await session.step(26, "Then \"Black\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(27, "And \"Asian\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(28, "When user presses ArrowUp", () => pressKey(page, "ArrowUp"));
      await session.step(29, "Then \"Asian\" legend item in legend of scatter plot viewer should be selected", () => shouldBe(page, el("\"Asian\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(30, "And \"Black\" legend item in legend of scatter plot viewer should not be selected", () => shouldNotBe(page, el("\"Black\" legend item in legend of scatter plot viewer"), "selected"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    }, 26);
  });
});
