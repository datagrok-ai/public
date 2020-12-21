class SDTMPackage extends DG.Package {
  //tags: app
  startApp(context) {
    let parser = document.createElement('a');
    parser.href = window.location;
    let pathSegments = parser.pathname.split('/');
    let urlLayout = (pathSegments.length > 3) ? decodeURI(pathSegments[3]) : null;

    let emptyTable = DG.DataFrame.create();

    let loadingLayout = false;
    let preview = ui.boolInput('', true);
    let limit = ui.intInput('Limit', 1000);
    let analysis = ui.div([], 'pure-form pure-form-aligned');
    let filtersHost = ui.div([analysis, ui.divH([preview.root, limit.root], 'sdtm-limit')], 'sdtm-controls pure-form');
    let count = ui.divH([ui.divText('Count:'), ui.loader()], 'sdtm-count');
    filtersHost.appendChild(count);

    let view = grok.shell.addTableView(emptyTable);
    view.name = 'SDTM:LB';
    view.basePath = '';
    view.description = 'SDTM LB domain viewer';
    view.root.className = 'grok-view grok-table-view sdtm-result';

    function updatePreview(t) {
      t.name = 'stdmlb';
      view.dataFrame = t;
    }

    updatePreview(emptyTable);
    ui.setUpdateIndicator(view.root, true);

    function selector(query, queryParameters, itemsColumn, caption, values, updateHandler) {
      let host = ui.divV([], 'sdtm-selector-host');
      host.appendChild(ui.loader());
      grok.data.query(query, queryParameters, true, 100).then(t => {
        let items = t.getCol(itemsColumn).toList();
        let input = ui.columnsInput(caption, DG.DataFrame.fromColumns(items.map(item => DG.Column.string(item))));
        if (values.length !== 0) {
          let _values = [];
          for (let value of values)
            if (items.includes(value))
              _values.push(value);
          input.stringValue = _values.join(',');
        }
        SDTMPackage.removeChildren(host);
        host.appendChild(input.root);

        input.onChanged(() => {
          updateHandler(input.stringValue.split(','));
        });

        updateHandler(input.stringValue.split(','));
      });
      return host;
    }

    let studiesParam = [];
    let testsParam = [];
    let sexParam = ['M', 'F'];
    let raceParam = ['OTHER', 'SPANISH, HISPANIC, OR LATINO', 'WHITE', 'MULTIPLE', 'BLACK, AFRICAN AMERICAN, OR NEGRO'];
    let methodParam = [];

    function aggregate(t) {
      grok.shell.setVar('table', t);
      return grok.functions.scriptSync('Aggregate(table, ' +
        'pivots = ["test"], ' +
        'aggregations = ["avg(lbstresn)"], ' +
        'groupByFields = ["studyid", "usubjid", "lbdy", "method", "units", "race", "sex"])');
    }

    function addQuotes(str) {
      return str.includes(' ') ? `"${str}"` : str;
    }

    function paramToInPattern(param) {
      return (param.length > 0) ? 'in (' + param.map(addQuotes).join(',') + ')' : null;
    }

    function queryUpdateResult(handler = () => {
    }) {
      if (loadingLayout)
        return;

      ui.setUpdateIndicator(view.root, true);

      let params = {
        'studies': paramToInPattern(studiesParam),
        'tests': paramToInPattern(testsParam),
        'sex': paramToInPattern(sexParam),
        'race': paramToInPattern(raceParam),
        'method': paramToInPattern(methodParam),
        'limit': preview.value ? limit.value : null
      };

      grok.data.query('SDTM:Results', params, true, 100).then(t => {
        updatePreview(aggregate(t));
        ui.setUpdateIndicator(view.root, false);
        handler();
      });

      SDTMPackage.removeChildren(count);
      count.appendChild(ui.divText('Count:'));
      count.appendChild(ui.loader());
      grok.data.query('SDTM:ResultsCount', params, true, 100).then(t => {
        SDTMPackage.removeChildren(count);
        count.appendChild(ui.divText(`Count: ${t.col('count').get(0)}`));
      });
    }

    let sex = ui.multiChoiceInput('Sex', sexParam, sexParam);
    sex.onChanged(() => {
      if (loadingLayout) return;
      sexParam = sex.value;
      queryUpdateResult();
    });
    filtersHost.appendChild(sex.root);

    let race = ui.multiChoiceInput('Race', raceParam, raceParam);
    race.onChanged(() => {
      if (loadingLayout) return;
      raceParam = race.value;
      queryUpdateResult();
    });
    filtersHost.appendChild(race.root);

    function removeFilters(after) {
      let idx = Array.from(filtersHost.childNodes).indexOf(after);
      if (idx !== -1) {
        while (filtersHost.childNodes.length > idx + 1)
          filtersHost.removeChild(filtersHost.childNodes[idx + 1]);
      }
    }

    function addFilters(updateResults = true) {
      removeFilters(race.root);
      let studiesSelector = selector('SDTM:Studies', {
        'race': paramToInPattern(raceParam)
      }, 'study', 'Studies', studiesParam, (studies) => {
        studiesParam = studies;
        removeFilters(studiesSelector);
        let testsSelector = selector('SDTM:Tests', {
          'studies': paramToInPattern(studiesParam),
          'race': paramToInPattern(raceParam)
        }, 'test', 'Tests', testsParam, (tests) => {
          testsParam = tests;
          removeFilters(testsSelector);
          let methodSelector = selector('SDTM:Methods', {
            'studies': paramToInPattern(studiesParam),
            'tests': paramToInPattern(testsParam),
            'race': paramToInPattern(raceParam)
          }, 'method', 'Methods', methodParam, (method) => {
            methodParam = method;
            if (updateResults)
              queryUpdateResult();
          });
          filtersHost.appendChild(methodSelector);
        });
        filtersHost.appendChild(testsSelector);
      });
      filtersHost.appendChild(studiesSelector);
    }

    preview.onChanged(() => {
      limit.enabled = preview.value;
      queryUpdateResult();
    });

    limit.onChanged(() => queryUpdateResult());


    // Analysis
    const STORAGE_NAME = 'sdtm-demo';
    let currentUserStorage = false;

    let nameInput = ui.stringInput('', 'Default');
    nameInput.captionLabel.style.width = '0';

    function saveToStorage() {
      let parameters = {
        'layout': view.saveLayout().toJson(),
        'studies': studiesParam,
        'tests': testsParam,
        'sex': sexParam,
        'race': raceParam,
        'method': methodParam
      };
      let layout = nameInput.value;
      grok.dapi.userDataStorage.postValue(STORAGE_NAME, layout, JSON.stringify(parameters), currentUserStorage);
      grok.shell.info(`Analysis "${layout}" is saved`);
      view.path = `/${layout}`;
    }

    function loadFromStorage(name, parameters) {
      view.path = `/${name}`;

      loadingLayout = true;

      parameters = JSON.parse(parameters);
      sexParam = parameters['sex'];
      sex.value = sexParam;
      raceParam = parameters['race'];
      race.value = raceParam;
      studiesParam = parameters['studies'];
      testsParam = parameters['tests'];
      methodParam = parameters['method'];

      addFilters(false);
      loadingLayout = false;
      queryUpdateResult(() => {
        view.loadLayout(DG.ViewLayout.fromJson(parameters['layout']));
      });
    }

    let save = ui.button('Save', saveToStorage);

    let load = ui.iconFA('bars', () => {
      grok.dapi.userDataStorage.get(STORAGE_NAME, currentUserStorage).then((layouts) => {
        if (layouts !== null && Object.keys(layouts).length === 0)
          grok.shell.info('Storage is empty. Save some analysis to the storage');
        else {
          let menu = DG.Menu.popup();
          for (let layout of Object.keys(layouts)) {
            menu.item(layout, () => {
              loadFromStorage(layout, layouts[layout]);
              nameInput.stringValue = layout;
            });
          }
          menu.show();
        }
      });
    }, 'Choose analysis from storage');

    let remove = ui.iconFA('trash-alt', () => {
      if (nameInput.value !== '') {
        grok.dapi.userDataStorage.remove(STORAGE_NAME, nameInput.value, currentUserStorage);
        grok.shell.info(`Analysis ${nameInput.value} is removed`);
        nameInput.value = 'Default';
      }
    }, 'Remove analysis from storage');

    /*
    let clear = ui.iconFA('trash-alt', () => {
        grok.dapi.userDataStorage.remove(STORAGE_NAME, null, currentUserStorage);
        grok.shell.info('Storage is cleared');
    }, 'Clear analysis storage');
    */

    let saveDiv = ui.divH([save, remove, /*clear*/]);
    saveDiv.style.marginLeft = '12px';

    analysis.appendChild(ui.divH([nameInput.root, load]));
    analysis.appendChild(saveDiv);

    let acc = view.toolboxPage.accordion;
    acc.addPane('SDTM', () => filtersHost, true, acc.panes[0]);

    grok.dapi.userDataStorage.get(STORAGE_NAME, currentUserStorage).then((layouts) => {
      if (layouts !== null && Object.keys(layouts).length === 0)
        addFilters();
      else {
        let layout = (urlLayout !== null && layouts.hasOwnProperty(urlLayout))
          ? urlLayout : Object.keys(layouts)[0];
        loadFromStorage(layout, layouts[layout]);
        nameInput.stringValue = layout;
      }
    });
  }

  //name: Medical History
  //description: Medical History info panel
  //tags: panel, widgets
  //input: string study {semType: Study}
  //output: widget result
  medicalHistoryStudy(study) {
    return this.getMedicalHistory(study, grok.shell.v.dataFrame.currentRow.usubjid);
  }

  //name: Medical History
  //description: Medical History info panel
  //tags: panel, widgets
  //input: string subj {semType: Subject}
  //output: widget result
  medicalHistorySubject(subj) {
    return this.getMedicalHistory(grok.shell.v.dataFrame.currentRow.studyid, subj);
  }

  //description: Gets medical history widget
  getMedicalHistory(study, subj) {
    let host = ui.div([ui.loader()]);
    grok.data.query('SDTM:MedicalHistory', {
      'study': study,
      'subj': subj
    }, true, 100).then(mi => {
      let grid = DG.Grid.create(mi);
      grid.root.style.width = '350px';
      grid.root.style.height = '400px';
      SDTMPackage.removeChildren(host);
      host.appendChild(grid.root);
    });
    return new DG.Widget(host);
  }

  //description: Removes all children from node
  static removeChildren(node) {
    while (node.firstChild)
      node.removeChild(node.firstChild);
  }
}
