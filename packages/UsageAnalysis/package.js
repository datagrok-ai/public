class UsageAnalysisPackage extends DG.Package {

  getSummary(summaryDf, counterType, title) {
    summaryDf.rows.match({counter_type: counterType}).select();
    let host = ui.block([],'d4-item-card card');
    host.appendChild(ui.h1(title));
    let list = [
      ['Day', summaryDf.get('day', summaryDf.selection.findNext(-1, 1))],
      ['Week', summaryDf.get('week', summaryDf.selection.findNext(-1, 1))],
      ['Month', summaryDf.get('month', summaryDf.selection.findNext(-1, 1))]
    ];
    summaryDf.selection.setAll(false);
    host.appendChild(ui.div([ui.table(list, (item, idx) =>
            [`${item[0]}:`, item[1]]
        )]
    ));
    return host;
  }

  //name: Usage Analysis
  //tags: app
  async startApp() {
    let view = grok.shell.newView('Usage Analysis');
    view.root.style.overflow = 'auto';
    let results = ui.div(this.users.root);
    //container for users tag editor
    let filters = ui.divH([],'filters');

    let summaryDf = await grok.data.query('UsageAnalysis:NewUsersEventsErrors');

    // containers for widgets
    let usersSummaryCard = this.getSummary(summaryDf, 'users_count', 'New users');
    let eventsSummaryCard = this.getSummary(summaryDf, 'events_count', 'New events');
    let errorsSummaryCard = this.getSummary(summaryDf, 'errors_count', 'New errors');
    let usageCard = null;
    let errorListCard = null;
    let uniqueUsersCard = null;
    let uniqueUsersPerDayCard = null;
    let eventTypeCard = null;
    let errorTypeCard = null;
    let testTrackingCard = null;

    let accToolbox = ui.accordion();
    let newUsersToolbox = ui.boolInput('New users', true, ()=>updateLayout());
    let newEventsToolbox = ui.boolInput('New events', true, ()=>updateLayout());
    let newErrorsToolbox = ui.boolInput('New errors', true, ()=>updateLayout());
    let uniqueUsersToolbox = ui.boolInput('Unique users', true, ()=>updateLayout());
    let uniqueUsersperDayToolbox = ui.boolInput('Unique users per day', true, ()=>updateLayout());
    uniqueUsersperDayToolbox.setTooltip('Unique users per day');
    let usageToolbox = ui.boolInput('Usage', true, ()=>updateLayout());
    let eventTypeToolbox = ui.boolInput('Event type', true, ()=>updateLayout());
    let errorTypeToolbox = ui.boolInput('Error type', true, ()=>updateLayout());
    let testTrackingToolbox = ui.boolInput('Test tracking', true, ()=>updateLayout());

    let date = ui.stringInput('Date', 'today');
    date.addPatternMenu('datetime');

    let userlist = null;
    let dock1 = grok.shell.dockManager;

    let addUser = ui.iconFA('plus', (x)=> {
      dock1.close(userlist);
      dock1.dock(userlist, 'left', null, 'Users', 0.17)
    });

    let users = this.users;
    let events = this.events;
    let isExactly = this.isExactly;

    async function showUsage() {

      uniqueUsersCard = ui.block([],'cardbox');
      uniqueUsersPerDayCard = ui.block([],'cardbox');
      usageCard = ui.block([],'cardbox');
      errorListCard = ui.block([],'cardbox');
      eventTypeCard = ui.block([],'cardbox');
      errorTypeCard = ui.block([],'cardbox');
      testTrackingCard = ui.block([],'cardbox');

      function addCardUsingDataframe(cardName, dataframe, f, supportUsers = true) {
        let host = ui.block([],'d4-item-card card');
        host.appendChild(ui.h1(cardName));

        if (cardName === 'Errors')
          grok.data.detectSemanticTypes(dataframe);
        host.appendChild(f(dataframe));
        return host;
      }

      function getCurrentFilter() {
        let filter = {'date': date.value, 'isExactly': isExactly.value};
        let selectedUsers = users.tags;

        filter['users'] = (selectedUsers.length !== 0) ? selectedUsers : ['all'];
        filter['events'] = (events.value.length !== 0) ? [isExactly.value ? events.value : events.value + '%'] : ['all'];

        return filter;
      }

      function addCardWithFilters(cardName, queryName, f) {
        let host = ui.block([],'d4-item-card card');
        host.appendChild(ui.h1(cardName));
        let loader = ui.loader();
        host.appendChild(loader);

        let filter = getCurrentFilter();

        grok.data.query('UsageAnalysis:' + queryName, filter).then((t) => {
          if (cardName === 'Errors')
            grok.data.detectSemanticTypes(t);
          host.appendChild(f(t));
          host.removeChild(loader);
        });
        return host;
      }

      while (results.firstChild)
        results.removeChild(results.firstChild);

      function subscribeOnTableWithEvents(table) {
        table.onCurrentRowChanged.subscribe((_) => {
          grok.dapi.log.include('session.user').find(table.currentRow.event_id).then((event) => grok.shell.o = event);
        });
      }

      function getUniqueUsers(cardName){
        let host = ui.block([],'d4-item-card card');
        host.style.overflow="auto";
        host.style.height="94.5%";
        host.appendChild(ui.h1(cardName));
        let loader = ui.loader();
        host.appendChild(loader);
        grok.data.query('UsageAnalysis:UniqueUsersByDate', {'date': date.value})
          .then(t => {
            let ids = Array.from(t.getCol('id').values());
            grok.dapi.getEntities(ids).then((users) => {
              host.appendChild(ui.list(users));
            });
          });
        host.removeChild(loader);
        return host;
      }

      let filter = getCurrentFilter();
      let errors = await grok.data.query('UsageAnalysis:ErrorsOnDateAndUsersAndEvents', filter);

      uniqueUsersCard.appendChild(getUniqueUsers('Unique Users'));
      uniqueUsersPerDayCard.appendChild(addCardWithFilters('Unique Users Per Day', 'UniqueUsersPerDayOnDateAndUsersAndEvents', (t) => DG.Viewer.lineChart(t).root));
      errorListCard.appendChild(addCardUsingDataframe('Error Types List', errors, (t) => DG.Viewer.grid(t).root));
      usageCard.appendChild(
            addCardWithFilters('Usage', 'EventsOnDateAndUsersAndEvents', (t) => {
              subscribeOnTableWithEvents(t);
              return DG.Viewer.scatterPlot(t, {'color': 'user'}).root;
            })
      );
      eventTypeCard.appendChild(addCardWithFilters('Event Types', 'EventsSummaryOnDate', (t) => DG.Viewer.barChart(t, {valueAggrType: 'avg'}).root));
      errorTypeCard.appendChild(addCardUsingDataframe('Error Types', errors, (t) => DG.Viewer.barChart(t, {valueAggrType: 'avg'}).root));
      testTrackingCard.appendChild(addCardWithFilters('Test Tracking', 'ManualActivityByDate', (t) => DG.Viewer.grid(t).root, false));
      results.append(
          ui.divH([usersSummaryCard,eventsSummaryCard,errorsSummaryCard]),
          ui.divH([uniqueUsersCard,uniqueUsersPerDayCard]),
          ui.divH([usageCard]),
          ui.divH([eventTypeCard,errorTypeCard]),
          ui.divH([errorListCard]),
          ui.divH([testTrackingCard])
      );

      updateUsersList();
      updateLayout();
    }

    date.onChanged(this.debounce(showUsage, 750));
    this.users.onChanged(showUsage);
    this.events.onChanged(this.debounce(showUsage, 750));
    this.isExactly.onChanged(showUsage);

    await showUsage();

    accToolbox.addPane('Services', () => ui.wait(async () => {
      let root = ui.div();
      let serviceInfos = await grok.dapi.admin.getServiceInfos();
      root.appendChild(ui.table(serviceInfos, (item, idx) =>
        [`${item.key}:`, $(ui.span([item.status])).addClass(`grok-plugin-status-${item.status.toLowerCase()}`)[0]]));
      return root;
    }), true);

    accToolbox.addPane('Filters', () => ui.div([
      date.root,
      this.events.root,
      this.isExactly.root,
      ui.div([ui.divH([ui.label('Add users'),addUser])])
    ],'toolbox-inputs'), false);

    accToolbox.addPane('Widgets', () => ui.inputs([
      newUsersToolbox,
      newEventsToolbox,
      newErrorsToolbox,
      uniqueUsersToolbox,
      uniqueUsersperDayToolbox,
      usageToolbox,
      eventTypeToolbox,
      errorTypeToolbox,
      testTrackingToolbox,
    ]), false);

    // Checking if the layout widget was hide or show
    function updateLayout(){
      (!newUsersToolbox.value ? $(usersSummaryCard).hide() : $(usersSummaryCard).show());
      (!newEventsToolbox.value ? $(eventsSummaryCard).hide() : $(eventsSummaryCard).show());
      (!newErrorsToolbox.value ? $(errorsSummaryCard).hide() : $(errorsSummaryCard).show());
      (!uniqueUsersToolbox.value ? $(uniqueUsersCard).hide() : $(uniqueUsersCard).show());
      (!uniqueUsersperDayToolbox.value ? $(uniqueUsersPerDayCard).hide() : $(uniqueUsersPerDayCard).show());
      (!usageToolbox.value ? $(usageCard).hide() : $(usageCard).show());
      (!eventTypeToolbox.value ? $(eventTypeCard).hide() : $(eventTypeCard).show());
      (!errorTypeToolbox.value ? $(errorTypeCard).hide() : $(errorTypeCard).show());
      (!testTrackingToolbox.value ? $(testTrackingCard).hide() : $(testTrackingCard).show());
    }

    // Add users to dock manager
    function updateUsersList(){
    userlist = ui.panel();
    $(userlist).addClass('user-list');
    grok.dapi.users
	   .order('login')
	    .list()
	    .then((allUsers) => userlist.appendChild(ui.divV(allUsers.map((p) => addUsersToList(p.name, p.login, p.picture)  ))));
    }

    // Checking if the tag editor doesn't have users and show them
    function addUsersToList(name,user,avatar){
      let className = '.'+name;
      let host = ui.divH([avatar, name, ui.button(ui.iconFA('plus-circle',()=>addUserToFilter(name, user))) ],name);
      $(host).css('align-items','center');
      if (users.tags.includes(user) == false){
        $(className).show()
        return host;
      }
    }
    // Add user from dock to tag editor and hide them from dock manager
    function addUserToFilter(name, user) {
      let className = '.'+name;
      users.addTag(user);
      return $(className).hide()
    }

    view.root.append(this.users.root);
    view.root.append(results);
    view.toolbox = accToolbox.root;
  }

  async init() {
    this.users = DG.TagEditor.create();
    this.events = ui.stringInput('Events', '');
    this.isExactly = ui.boolInput('Exact m...', false);
    this.isExactly.setTooltip('Exact matching');

    grok.events.onContextMenu.subscribe((args) => {
      if (args.args.item instanceof DG.User) {
        args.args.menu.item('Add user to filter', async () => this.addUserToFilter(args.args.item.login));
      } else if (args.args.item instanceof DG.LogEvent) {
        args.args.menu.item('Add event to filter', async () => {
          this.events.value = args.args.item.name;
          this.isExactly.value = true;
        });
      }
    });
  }

  debounce(fn, time) {
    let timeout;
    return function () {
      const functionCall = () => fn.apply(this, arguments);
      clearTimeout(timeout);
      timeout = setTimeout(functionCall, time);
    }
  }


  //name: Create JIRA ticket
  //description: Creates JIRA ticket using current error log
  //tags: panel, widgets
  //input: string msg {semType: ErrorMessage}
  //output: widget result
  //condition: true
  createJiraTicket(msg) {
    let root = ui.div();

    let summary = ui.stringInput('Summary', '');
    let description = ui.stringInput('Description', msg);

    let button = ui.bigButton('CREATE', () => {
      grok.data.query('Vnerozin:JiraCreateIssue', {
        'createRequest': JSON.stringify({
          "fields": {
            "project": {
              "key": "GROK"
            },
            "summary": summary.value,
            "description": description.value,
            "issuetype": {
              "name": "Bug"
            }
          }
        }),
        'updateHistory': false,
      }).then((t) => {
        grok.shell.info('Created');
      });
    });
    button.style.marginTop = '12px';

    root.appendChild(ui.inputs([summary, description]));
    root.appendChild(button);

    return new DG.Widget(root);
  }
}
