import "dart:async";
import 'dart:io';

import 'package:args/args.dart';
import 'package:path/path.dart' as path;

Future main(List<String> args) async {
  var parser = new ArgParser()
    ..addOption('dest_path',
        abbr: 'd',
        defaultsTo: '.',
        help: 'Path to directory to save generated file')
    ..addOption('source_path',
        abbr: 's',
        defaultsTo: '.',
        help: 'Path to directory from which generate the file')
    ..addFlag('help', abbr: 'h', help: 'Show help');

  var params = parser.parse(args);

  if (params['help']) {
    print(parser.usage);
    return;
  }

  Future<Map> getTree(String rootPath) async {
    var filesPaths = await new Directory(rootPath)
        .list(recursive: true)
        .where((fse) => (fse is File && path.extension(fse.path) == '.html'))
        .map((m) => path.posix.joinAll(path.split(m.path)))
        .toList();

    Map map = {};
    for (int n = 0; n < filesPaths.length; n++) {
      var filePath = filesPaths[n];
      var relativePath = path.posix.relative(filePath, from: rootPath);
      var content = await (new File(filePath)).readAsString();
      String title = new RegExp(r'<!--( ?)TITLE:( ?)(.*)-->')
              .firstMatch(content)
              ?.group(3)
              ?.trim() ??
          path.split(relativePath).last.replaceAll('.html', '');
      var _map = map;
      for (var p in path.split(relativePath)) {
        if (p.endsWith('.html'))
          _map[title] = {'url': '/help/' + relativePath, 'title': title};
        else {
          if (!_map.containsKey(p)) _map[p] = {};
          _map = _map[p];
        }
      }
    }

    return map;
  }

  var links = await getTree(params['source_path']);

  String treeNodeHtml(Map links) {
    var treeHtml = '';
    for (var k in links.keys
        .toList()
        .where((k) => !links[k].containsKey('url'))
        .toList()
      ..sort()) {
      treeHtml +=
          '<li><label>${k[0].toUpperCase() + k.substring(1)}</label><ul>${treeNodeHtml(links[k])}</ul>';
    }
    for (var k in links.keys
        .toList()
        .where(
            (k) => links[k].containsKey('url') && !links[k]?.isEmpty ?? false)
        .toList()
      ..sort()) {
      treeHtml +=
          '<li><a href="${links[k]['url']}">${links[k]['title'][0].toUpperCase() + links[k]['title'].substring(1)}</a></li>';
    }
    return treeHtml;
  }

  var menuHtml = '<ul>${treeNodeHtml(links)}</ul>';
  File(params['dest_path'] + '/menu.html').writeAsString(menuHtml);
}
