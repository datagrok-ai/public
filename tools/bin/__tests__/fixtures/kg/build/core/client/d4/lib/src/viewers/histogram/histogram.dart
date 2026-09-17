import 'dart:html';
import 'package:d4/src/viewers/viewer.dart';
import '../legend_cache.dart' show LegendCache;
import 'package:quiver/iterables.dart';
import 'package:d4/src/viewers/nowhere.dart';

class Histogram extends Viewer {
  String get helpUrl => HelpUrl.Histogram;
}
