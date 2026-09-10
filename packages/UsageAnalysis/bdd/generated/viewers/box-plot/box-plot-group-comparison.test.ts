/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-group-comparison.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, selectIn, shouldBe, shouldContainText, shouldHaveText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, removeColumn, tableColumnComplete, tableColumns, tableOpen, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, hasArea, hasNoArea, hoverArea, pickFromAreaContextMenu, pointerAway, propertyShouldBe, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot group comparison", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-group-comparison.feature", import.meta.url);
  test("Box plot group comparison", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","SEX"]]));
    await run.scenario("The bare p-value and its reveal icon", async () => {
      await session.step(17, "Then \"Show P Value\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "true"));
      await session.step(18, "And \"Show Group Comparison\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "false"));
      await session.step(19, "And box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
      await session.step(20, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(21, "Then tooltip should contain text \"t-test\"", () => shouldContainText(page, el("tooltip"), "t-test"));
      await session.step(22, "And show-group-stats icon in box plot viewer should be visible", () => shouldBe(page, el("show-group-stats icon in box plot viewer"), "visible"));
      await session.step(23, "When user clicks on show-group-stats icon in box plot viewer", () => clickOn(page, el("show-group-stats icon in box plot viewer")));
      await session.step(24, "Then \"Show Group Comparison\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "true"));
      await session.step(25, "And box plot viewer should have a \"group comparison\" area", () => hasArea(page, el("box plot viewer"), "group comparison"));
    });
    await run.scenario("The overall test and the method selector", async () => {
      await session.step(28, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(29, "Then tooltip should contain text \"Welch\"", () => shouldContainText(page, el("tooltip"), "Welch"));
      await session.step(30, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(31, "And user selects \"Student\" in method choice input in box plot viewer", () => selectIn(page, "Student", el("method choice input in box plot viewer")));
      await session.step(32, "Then \"Method\" property of box plot viewer should be \"Student\"", () => propertyShouldBe(page, "Method", el("box plot viewer"), "Student"));
      await session.step(33, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(34, "Then tooltip should contain text \"Student\"", () => shouldContainText(page, el("tooltip"), "Student"));
      await session.step(35, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(36, "And user selects \"Welch\" in method choice input in box plot viewer", () => selectIn(page, "Welch", el("method choice input in box plot viewer")));
      await session.step(37, "Then \"Method\" property of box plot viewer should be \"Welch\"", () => propertyShouldBe(page, "Method", el("box plot viewer"), "Welch"));
    });
    await run.scenario("Three groups, a control group and its comparisons table", async () => {
      await session.step(40, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(41, "Then box plot viewer should have repainted", () => repainted(page, el("box plot viewer")));
      await session.step(42, "And box plot viewer should have a \"category Asian\" area", () => hasArea(page, el("box plot viewer"), "category Asian"));
      await session.step(43, "When user hovers over the \"p value\" area of box plot viewer", () => hoverArea(page, "p value", el("box plot viewer")));
      await session.step(44, "Then tooltip should contain text \"ANOVA\"", () => shouldContainText(page, el("tooltip"), "ANOVA"));
      await session.step(45, "When user moves the pointer away from box plot viewer", () => pointerAway(page, el("box plot viewer")));
      await session.step(46, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(47, "And user selects \"Caucasian\" in control-group choice input in box plot viewer", () => selectIn(page, "Caucasian", el("control-group choice input in box plot viewer")));
      await session.step(48, "Then \"Control Comparisons\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Control Comparisons", el("box plot viewer"), "true"));
      await session.step(49, "And \"Control Group\" property of box plot viewer should be \"Caucasian\"", () => propertyShouldBe(page, "Control Group", el("box plot viewer"), "Caucasian"));
      await session.step(50, "And box plot viewer should have a \"p value of Asian\" area", () => hasArea(page, el("box plot viewer"), "p value of Asian"));
      await session.step(51, "And box plot viewer should have a \"control band\" area", () => hasArea(page, el("box plot viewer"), "control band"));
      await session.step(52, "When user picks \"Add Control Comparisons Table\" from the context menu of the \"group comparison\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Add Control Comparisons Table", "group comparison", el("box plot viewer")));
      await session.step(53, "Then table \"Control Comparisons: AGE by RACE vs Caucasian\" should be open", () => tableOpen(page, "Control Comparisons: AGE by RACE vs Caucasian"));
      await session.step(54, "And table \"Control Comparisons: AGE by RACE vs Caucasian\" should have 3 rows", () => tableRows(page, "Control Comparisons: AGE by RACE vs Caucasian", 3));
      await session.step(55, "And table \"Control Comparisons: AGE by RACE vs Caucasian\" should have columns \"Conclusion, Group, n, Mean, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g\"", () => tableColumns(page, "Control Comparisons: AGE by RACE vs Caucasian", "Conclusion, Group, n, Mean, Mean diff, 95% CI low, 95% CI high, t, df, p (raw), p (adj), Hedges' g"));
      await session.step(56, "And table \"Control Comparisons: AGE by RACE vs Caucasian\" should have no missing values in \"p (adj)\" column", () => tableColumnComplete(page, "Control Comparisons: AGE by RACE vs Caucasian", "p (adj)"));
      await session.step(57, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
      await session.step(58, "And user clicks on the \"p value of Asian\" area of box plot viewer", () => clickArea(page, "p value of Asian", el("box plot viewer")));
      await session.step(59, "Then Results section in context panel should be present", () => shouldBe(page, el("Results section in context panel"), "present"));
      await session.step(60, "And Statistics section in context panel should contain text \"Asian\"", () => shouldContainText(page, el("Statistics section in context panel"), "Asian"));
    });
    await run.scenario("Two-way ANOVA", async () => {
      await session.step(63, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Control Comparisons","false"],["Control Group",""],["Category 2","SEX"]]));
      await session.step(67, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(68, "Then baseline choice input in box plot viewer should be visible", () => shouldBe(page, el("baseline choice input in box plot viewer"), "visible"));
      await session.step(69, "And box plot viewer should have a \"RACE effect\" area", () => hasArea(page, el("box plot viewer"), "RACE effect"));
      await session.step(70, "When user picks \"Add Two-Way ANOVA Table\" from the context menu of the \"group comparison\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Add Two-Way ANOVA Table", "group comparison", el("box plot viewer")));
      await session.step(71, "Then table \"Two-Way ANOVA: AGE by RACE, SEX\" should be open", () => tableOpen(page, "Two-Way ANOVA: AGE by RACE, SEX"));
      await session.step(72, "And table \"Two-Way ANOVA: AGE by RACE, SEX\" should have 5 rows", () => tableRows(page, "Two-Way ANOVA: AGE by RACE, SEX", 5));
      await session.step(73, "And table \"Two-Way ANOVA: AGE by RACE, SEX\" should have columns \"Conclusion, Source of variance, SS, DF, MS, F, p-value\"", () => tableColumns(page, "Two-Way ANOVA: AGE by RACE, SEX", "Conclusion, Source of variance, SS, DF, MS, F, p-value"));
      await session.step(74, "And table \"Two-Way ANOVA: AGE by RACE, SEX\" should have no missing values in \"SS\" column", () => tableColumnComplete(page, "Two-Way ANOVA: AGE by RACE, SEX", "SS"));
      await session.step(75, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    });
    await run.scenario("Closing the comparison", async () => {
      await session.step(78, "When user hovers over the \"group comparison\" area of box plot viewer", () => hoverArea(page, "group comparison", el("box plot viewer")));
      await session.step(79, "And user clicks on close-group-stats icon in box plot viewer", () => clickOn(page, el("close-group-stats icon in box plot viewer")));
      await session.step(80, "Then \"Show Group Comparison\" property of box plot viewer should be \"false\"", () => propertyShouldBe(page, "Show Group Comparison", el("box plot viewer"), "false"));
      await session.step(81, "And \"Show P Value\" property of box plot viewer should be \"true\"", () => propertyShouldBe(page, "Show P Value", el("box plot viewer"), "true"));
      await session.step(82, "And method choice input in box plot viewer should be hidden", () => shouldBe(page, el("method choice input in box plot viewer"), "hidden"));
      await session.step(83, "And box plot viewer should not have a \"group comparison\" area", () => hasNoArea(page, el("box plot viewer"), "group comparison"));
      await session.step(84, "And box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
    });
    await run.scenario("A covariate adjusts the value axis", async () => {
      await session.step(87, "Given user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 2",""],["Control Comparisons","false"],["Control Group",""],["Category 1","SEX"],["Value","WEIGHT"],["Show Group Comparison","true"]]));
      await session.step(94, "When user sets \"Adjust By\" property of box plot viewer to \"HEIGHT\"", () => setProperty(page, "Adjust By", el("box plot viewer"), "HEIGHT"));
      await session.step(95, "Then \"Adjustment\" property of box plot viewer should be \"regressOut\"", () => propertyShouldBe(page, "Adjustment", el("box plot viewer"), "regressOut"));
      await session.step(96, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(97, "Then \"Adjust by\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Adjust by\" column input in box plot viewer"), "visible"));
      await session.step(98, "And \"Adjust by\" column input in box plot viewer should contain text \"HEIGHT\"", () => shouldContainText(page, el("\"Adjust by\" column input in box plot viewer"), "HEIGHT"));
      await session.step(99, "And Value column input in box plot viewer should have text \"WEIGHT\"", () => shouldHaveText(page, el("Value column input in box plot viewer"), "WEIGHT"));
      await session.step(100, "And adjustment choice input in box plot viewer should have the value \"Regress-out\"", () => shouldHaveValue(page, el("adjustment choice input in box plot viewer"), "Regress-out"));
      await session.step(101, "When user selects \"Ratio\" in adjustment choice input in box plot viewer", () => selectIn(page, "Ratio", el("adjustment choice input in box plot viewer")));
      await session.step(102, "Then \"Adjustment\" property of box plot viewer should be \"ratio\"", () => propertyShouldBe(page, "Adjustment", el("box plot viewer"), "ratio"));
    });
    await run.scenario("ANCOVA against a control group", async () => {
      await session.step(105, "When user sets \"Category 1\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Category 1", el("box plot viewer"), "RACE"));
      await session.step(106, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(107, "And user selects \"Caucasian\" in control-group choice input in box plot viewer", () => selectIn(page, "Caucasian", el("control-group choice input in box plot viewer")));
      await session.step(108, "And user selects \"ANCOVA\" in method choice input in box plot viewer", () => selectIn(page, "ANCOVA", el("method choice input in box plot viewer")));
      await session.step(109, "Then \"Method\" property of box plot viewer should be \"ANCOVA\"", () => propertyShouldBe(page, "Method", el("box plot viewer"), "ANCOVA"));
      await session.step(110, "And \"Adjustment\" property of box plot viewer should be \"ratio\"", () => propertyShouldBe(page, "Adjustment", el("box plot viewer"), "ratio"));
      await session.step(111, "And adjustment choice input in box plot viewer should be hidden", () => shouldBe(page, el("adjustment choice input in box plot viewer"), "hidden"));
      await session.step(112, "When user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(113, "Then \"Adjust by\" column input in box plot viewer should be visible", () => shouldBe(page, el("\"Adjust by\" column input in box plot viewer"), "visible"));
      await session.step(114, "And \"Adjust by\" column input in box plot viewer should contain text \"HEIGHT\"", () => shouldContainText(page, el("\"Adjust by\" column input in box plot viewer"), "HEIGHT"));
      await session.step(115, "When user picks \"Add ANCOVA Table\" from the context menu of the \"group comparison\" area of box plot viewer", () => pickFromAreaContextMenu(page, "Add ANCOVA Table", "group comparison", el("box plot viewer")));
      await session.step(116, "Then table \"ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian\" should be open", () => tableOpen(page, "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian"));
      await session.step(117, "And table \"ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian\" should have 4 rows", () => tableRows(page, "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian", 4));
      await session.step(118, "And table \"ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian\" should have columns \"Conclusion, Group, n, Raw mean, Adjusted mean, SE, Adj. diff, p-value, Hedges' g\"", () => tableColumns(page, "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian", "Conclusion, Group, n, Raw mean, Adjusted mean, SE, Adj. diff, p-value, Hedges' g"));
      await session.step(119, "And table \"ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian\" should have no missing values in \"Adjusted mean\" column", () => tableColumnComplete(page, "ANCOVA: WEIGHT (ANCOVA adj) by RACE vs Caucasian", "Adjusted mean"));
      await session.step(120, "When user switches to the \"demog-1000\" table view", () => switchTableView(page, "demog-1000"));
    });
    await run.scenario("The matched baseline and the Simpson's paradox cue", async () => {
      await session.step(123, "When user sets \"Category 2\" property of box plot viewer to \"SEX\"", () => setProperty(page, "Category 2", el("box plot viewer"), "SEX"));
      await session.step(124, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(125, "And user selects \"Matched · per stratum\" in baseline choice input in box plot viewer", () => selectIn(page, "Matched · per stratum", el("baseline choice input in box plot viewer")));
      await session.step(126, "Then \"Baseline Mode\" property of box plot viewer should be \"matched\"", () => propertyShouldBe(page, "Baseline Mode", el("box plot viewer"), "matched"));
      await session.step(127, "And box plot viewer should have a \"control band F\" area", () => hasArea(page, el("box plot viewer"), "control band F"));
      await session.step(128, "And box plot viewer should have a \"control band M\" area", () => hasArea(page, el("box plot viewer"), "control band M"));
      await session.step(129, "And box plot viewer should not have a \"control band\" area", () => hasNoArea(page, el("box plot viewer"), "control band"));
      await session.step(130, "When user adds a calculated column \"SIMPSON_STRAT\" with formula \"if(Mod(Round(${HEIGHT} * 1000), 2) == 0, \\\"A\\\", \\\"B\\\")\"", () => addCalculated(page, "SIMPSON_STRAT", "if(Mod(Round(${HEIGHT} * 1000), 2) == 0, \"A\", \"B\")"));
      await session.step(131, "And user adds a calculated column \"SIMPSON_VAL\" with formula \"if(${SEX} == \\\"M\\\", 0, if(Mod(Round(${HEIGHT} * 1000), 2) == 0, 2.5, -2.5)) + (Mod(Round(${WEIGHT} * 137), 600) / 30)\"", () => addCalculated(page, "SIMPSON_VAL", "if(${SEX} == \"M\", 0, if(Mod(Round(${HEIGHT} * 1000), 2) == 0, 2.5, -2.5)) + (Mod(Round(${WEIGHT} * 137), 600) / 30)"));
      await session.step(132, "And user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1","SEX"],["Category 2","SIMPSON_STRAT"],["Value","SIMPSON_VAL"]]));
      await session.step(136, "And user hovers over box plot viewer", () => hoverOver(page, el("box plot viewer")));
      await session.step(137, "And user selects \"M\" in control-group choice input in box plot viewer", () => selectIn(page, "M", el("control-group choice input in box plot viewer")));
      await session.step(138, "Then box plot viewer should have a \"control band A\" area", () => hasArea(page, el("box plot viewer"), "control band A"));
      await session.step(139, "And simpson-warning icon in box plot viewer should be visible", () => shouldBe(page, el("simpson-warning icon in box plot viewer"), "visible"));
      await session.step(140, "When user hovers over simpson-warning icon in box plot viewer", () => hoverOver(page, el("simpson-warning icon in box plot viewer")));
      await session.step(141, "Then tooltip should contain text \"Pooling cancels opposite within-stratum trends\"", () => shouldContainText(page, el("tooltip"), "Pooling cancels opposite within-stratum trends"));
      await session.step(142, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Value","WEIGHT"],["Category 1","RACE"],["Category 2","SEX"]]));
      await session.step(146, "And user removes \"SIMPSON_STRAT\" column", () => removeColumn(page, "SIMPSON_STRAT"));
      await session.step(147, "And user removes \"SIMPSON_VAL\" column", () => removeColumn(page, "SIMPSON_VAL"));
      await session.step(148, "And user sets \"Adjust By\" property of box plot viewer to \"\"", () => setProperty(page, "Adjust By", el("box plot viewer"), ""));
      await session.step(149, "Then \"Adjust By\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Adjust By", el("box plot viewer"), ""));
      await session.step(150, "And \"Adjustment\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Adjustment", el("box plot viewer"), ""));
      await session.step(151, "And \"Method\" property of box plot viewer should be \"\"", () => propertyShouldBe(page, "Method", el("box plot viewer"), ""));
    });
    run.finish();
  });
});
