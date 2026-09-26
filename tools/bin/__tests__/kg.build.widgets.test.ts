/// Heritage from the lexical Dart pass and the widget nodes it yields (d4-features.md §4 item A): `extends` and
/// `implements` between the types of a package and the ones it imports, an ambiguous base drawn for nothing, and
/// one `widget` node per class whose chain reaches d4's Widget, wherever the class sits.
import {describe, it, expect} from 'vitest';
import {copyFixture, buildFixture, write} from './kg-fixture';

const WIDGETS = 'core/client/d4/lib/src/widgets';
const WIDGET = `${WIDGETS}/widget.dart`;
const INPUT_BASE = `${WIDGETS}/inputs/input_base.dart`;
const TEXT_INPUT = `${WIDGETS}/inputs/text_input.dart`;
const DIALOG = `${WIDGETS}/dialog/dialog.dart`;
const BALLOON = `${WIDGETS}/balloon/balloon.dart`;
const MENU_POPUP = `${WIDGETS}/menu/popup.dart`;
const TOOLTIP_POPUP = `${WIDGETS}/tooltip/popup.dart`;
const FILTER_CONTROL = 'core/client/d4/lib/src/viewers/filters/filter_control.dart';
const CONSOLE = 'core/client/xamgle/lib/src/features/console.dart';
const VIEW_BASE = `${WIDGETS}/view_base.dart`;
const VIEWER_BASE = 'core/client/d4/lib/src/viewer_base/viewer_base.dart';
const TABLE_VIEW = 'core/client/xamgle/lib/src/views/table_view.dart';
const PIE_CORE = 'core/client/d4/lib/src/viewers/pie_chart/pie_chart_core.dart';

const graph = buildFixture(copyFixture('build', (repo) => {
  write(repo, WIDGET, ['class PropMixin {}', 'mixin Detachable {}', 'abstract class IRoot {}', '', 'abstract class Widget extends PropMixin',
    '    with Detachable implements IRoot, IFunctionsProvider {', '}', ''].join('\n'));
  write(repo, INPUT_BASE, ["import '../widget.dart';", '', 'abstract class InputBase<T> extends GrokJsObject {}', ''].join('\n'));
  write(repo, TEXT_INPUT, ["import 'input_base.dart';", '', 'class TextInput extends InputBase<String> {}', 'class TextFormatter {}', ''].join('\n'));
  write(repo, DIALOG, ["import 'package:d4/src/widgets/widget.dart';", "import 'dart:html';", '', 'class Dialog extends Widget {}', 'class HtmlWidget extends DivElement {}', ''].join('\n'));
  write(repo, MENU_POPUP, 'class Popup {}\n');
  write(repo, TOOLTIP_POPUP, "import '../widget.dart';\n\nclass Popup extends Widget {}\n");
  write(repo, BALLOON, "import '../widget.dart';\n\nclass Balloon extends Popup {}\n");
  write(repo, FILTER_CONTROL, "import '../../widgets/widget.dart';\n\nclass FilterControl extends Widget {}\n");
  write(repo, CONSOLE, "import 'package:d4/src/widgets/widget.dart';\n\nclass Console extends Widget {}\n");
  write(repo, VIEW_BASE, "import 'widget.dart';\n\nabstract class ViewBase extends Widget {}\n");
  write(repo, VIEWER_BASE, "import '../widgets/widget.dart';\n\nabstract class ViewerBase extends Widget {}\n");
  write(repo, TABLE_VIEW, "import 'package:d4/src/widgets/view_base.dart';\n\nclass TableView extends ViewBase {}\n");
  write(repo, PIE_CORE, "import '../../viewer_base/viewer_base.dart';\n\nclass PieChartCore extends ViewerBase {}\n");
}), 'homes,dart,membership');
const decl = (file: string, name: string) => `decl:${file}#${name}`;

describe('Dart heritage and widget nodes', () => {
  it('draws extends from a class to the one type its package or an imported one declares under the base name, type arguments stripped', async () => {
    const {rows, problems} = await graph;
    const extend = rows('edges/extends');
    expect(extend.map((e) => [e.from, e.to]).sort()).toEqual([
      [decl(CONSOLE, 'Console'), decl(WIDGET, 'Widget')],
      [decl(DIALOG, 'Dialog'), decl(WIDGET, 'Widget')],
      [decl(FILTER_CONTROL, 'FilterControl'), decl(WIDGET, 'Widget')],
      [decl(PIE_CORE, 'PieChartCore'), decl(VIEWER_BASE, 'ViewerBase')],
      [decl(TABLE_VIEW, 'TableView'), decl(VIEW_BASE, 'ViewBase')],
      [decl(TEXT_INPUT, 'TextInput'), decl(INPUT_BASE, 'InputBase')],
      [decl(TOOLTIP_POPUP, 'Popup'), decl(WIDGET, 'Widget')],
      [decl(VIEWER_BASE, 'ViewerBase'), decl(WIDGET, 'Widget')],
      [decl(VIEW_BASE, 'ViewBase'), decl(WIDGET, 'Widget')],
      [decl(WIDGET, 'Widget'), decl(WIDGET, 'PropMixin')],
      [decl('core/client/d4/lib/src/viewers/histogram/histogram.dart', 'Histogram'), decl('core/client/d4/lib/src/viewers/viewer.dart', 'Viewer')],
      [decl('core/client/d4/lib/src/viewers/scatterplot/scatter.dart', 'ScatterPlot'), decl('core/client/d4/lib/src/viewers/viewer.dart', 'Viewer')],
    ].sort());
    expect(extend[0]).toMatchObject({derived_by: 'lexical', confidence: 0.9, evidence: [extend[0].from.slice('decl:'.length).split('#')[0]]});
    expect(rows('edges/implements').map((e) => [e.from, e.to]).sort()).toEqual([
      [decl(WIDGET, 'Widget'), decl(WIDGET, 'Detachable')],
      [decl(WIDGET, 'Widget'), decl(WIDGET, 'IRoot')],
    ]);
    // a base two files declare is ambiguous; one nobody walked (dart:html) is nothing
    expect(problems.ambiguous_extends).toEqual([`${BALLOON}: Balloon extends Popup is declared in ${MENU_POPUP} and ${TOOLTIP_POPUP}`]);
    expect(extend.some((e) => e.from === decl(DIALOG, 'HtmlWidget'))).toBe(false);
  });

  // the roots themselves are not widgets, nor is anything below ViewBase (views) or ViewerBase (viewers)
  it('makes a widget node of every class whose extends chain reaches Widget or InputBase without passing a view or viewer base', async () => {
    const {rows} = await graph;
    const widgets = rows('nodes/widget');
    expect(widgets.map((w) => [w.id, w.path, w.base, w.abstract]).sort()).toEqual([
      ['widget:Console', CONSOLE, 'Widget', false],
      ['widget:Dialog', DIALOG, 'Widget', false],
      ['widget:FilterControl', FILTER_CONTROL, 'Widget', false],
      ['widget:Popup', TOOLTIP_POPUP, 'Widget', false],
      ['widget:TextInput', TEXT_INPUT, 'InputBase', false],
    ]);
    expect(widgets.find((w) => w.id === 'widget:TextInput')).toMatchObject({type: 'widget', name: 'TextInput', line: 3, declaration: decl(TEXT_INPUT, 'TextInput'),
      builtin: true, language: 'dart', provenance: 'lexical', source_layer: 'core'});
    expect(rows('edges/declaration').filter((e) => e.from.startsWith('widget:')).map((e) => [e.from, e.to]).sort()).toEqual(
      widgets.map((w) => [w.id, w.declaration]).sort());
    expect(rows('edges/declares').filter((e) => e.to.startsWith('widget:')).map((e) => [e.from, e.to]).sort()).toEqual(
      widgets.map((w) => [`file:${w.path}`, w.id]).sort());
    expect(rows('edges/declares').find((e) => e.to === 'widget:Dialog')).toMatchObject({derived_by: 'lexical', confidence: 0.9, evidence: [DIALOG]});
  });
});
