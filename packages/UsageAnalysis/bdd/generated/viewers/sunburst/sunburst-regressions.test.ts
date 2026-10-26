/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/sunburst/sunburst-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.sunburst, public#2292]
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
import {clickOn, selectIn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, noErrors, painted, readingAsRemembered, readingNotAsRemembered, rememberReading, reportsNoError} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Sunburst regressions of closed issues", () => {
  const session = feature(test, "features/viewers/sunburst/sunburst-regressions.feature", import.meta.url);
  test("The Table property moves the sunburst to another table and back (public#2292)", {tag: ["@viewers", "@realizes:charts.viewer.sunburst", "@realizes:public#2292"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a sunburst viewer", () => addViewer(page, "sunburst"));
    await session.step(18, "Then sunburst viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("sunburst viewer"), "demog-1000"));
    await session.step(19, "And sunburst viewer should be painted", () => painted(page, el("sunburst viewer")));
    await session.step(20, "When user remembers the \"segment names\" reading of sunburst viewer", () => rememberReading(page, "segment names", el("sunburst viewer")));
    await session.step(21, "And user clicks on settings icon of sunburst viewer", () => clickOn(page, el("settings icon of sunburst viewer")));
    await session.step(22, "And user selects \"spgi-100\" in \"Table\" property in context panel", () => selectIn(page, "spgi-100", el("\"Table\" property in context panel")));
    await session.step(23, "Then sunburst viewer should be bound to table \"spgi-100\"", () => boundTable(page, el("sunburst viewer"), "spgi-100"));
    await session.step(24, "And the \"segment names\" reading of sunburst viewer should not be as remembered", () => readingNotAsRemembered(page, "segment names", el("sunburst viewer")));
    await session.step(25, "And sunburst viewer should report no error", () => reportsNoError(page, el("sunburst viewer")));
    await session.step(26, "When user selects \"demog-1000\" in \"Table\" property in context panel", () => selectIn(page, "demog-1000", el("\"Table\" property in context panel")));
    await session.step(27, "Then sunburst viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("sunburst viewer"), "demog-1000"));
    await session.step(28, "And the \"segment names\" reading of sunburst viewer should be as remembered", () => readingAsRemembered(page, "segment names", el("sunburst viewer")));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
});
