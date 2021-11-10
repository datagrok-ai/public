import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import '../css/usage_analysis.css';

export class UsageAnalysisView extends DG.ViewBase {

    getSummary(summaryDf:any, counterType:string, title:string) {
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

    constructor(name: string) {
        super();
        this.name = name;

        // toolbox accordion 
        let accToolbox = ui.accordion();

        // toolbox inputs
        let dateInput = ui.stringInput('Date', 'today');
        dateInput.addPatternMenu('datetime');
        dateInput.setTooltip('Set the date period');

        let usersInput = ui.stringInput('User', 'all');
        usersInput.setTooltip('Enter users login');

        let eventsInput = ui.stringInput('Event', 'all');
        eventsInput.setTooltip('Enter event name');

        let isExactlyInput = ui.boolInput('Exact m...', true);
        isExactlyInput.setTooltip('Exact matching');

        let applyBtn = ui.bigButton('Apply',()=>{this.debounce(showUsage(), 750)});

        // containers for widgets
        let usersSummaryCard:HTMLDivElement;
        let eventsSummaryCard:HTMLDivElement;
        let errorsSummaryCard:HTMLDivElement;
        let usageCard = null;
        let errorListCard = null;
        let uniqueUsersCard = null;
        let uniqueUsersPerDayCard = null;
        let eventTypeCard = null;
        let errorTypeCard = null;
        let testTrackingCard = null;

        accToolbox.addPane('Services', () => ui.wait(async () => {
            let root = ui.div();
            let serviceInfos = await grok.dapi.admin.getServiceInfos();
            root.appendChild(ui.table(serviceInfos, (item, idx) =>
            //@ts-ignore
              [`${item.key}:`, $(ui.span([item.status])).addClass(`grok-plugin-status-${item.status.toLowerCase()}`)[0]]));
            return root;
          }), true);
      
        accToolbox.addPane('Filters', () => ui.div([
            dateInput,
            usersInput,
            eventsInput,
            isExactlyInput,
            applyBtn
        ],'toolbox-inputs'), true);

        let showUsage = async () => {
            ui.setUpdateIndicator(this.root, true);
            // summary cards
            let summaryDf = await grok.data.query('UsageAnalysis:NewUsersEventsErrors');
            usersSummaryCard = this.getSummary(summaryDf, 'users_count', 'New users');
            eventsSummaryCard = this.getSummary(summaryDf, 'events_count', 'New events');
            errorsSummaryCard = this.getSummary(summaryDf, 'errors_count', 'New errors');
            
            // usage carcs containers
            uniqueUsersCard = ui.block([],'cardbox');
            uniqueUsersPerDayCard = ui.block([],'cardbox');
            usageCard = ui.block([],'cardbox');
            errorListCard = ui.block([],'cardbox');
            eventTypeCard = ui.block([],'cardbox');
            errorTypeCard = ui.block([],'cardbox');
            testTrackingCard = ui.block([],'cardbox');
            
            while (this.root.firstChild)
                this.root.removeChild(this.root.firstChild);

            function subscribeOnTableWithEvents(table:any) {
                table.onCurrentRowChanged.subscribe((_: any) => {
                grok.dapi.log.include('session.user').find(table.currentRow.event_id).then((event) => grok.shell.o = event);
                });
            }
            
            let filter = getCurrentFilter();
            let errors = await grok.data.query('UsageAnalysis:ErrorsOnDateAndUsersAndEvents', filter);

            uniqueUsersCard.appendChild(getUniqueUsers('Unique Users'));
            uniqueUsersPerDayCard.appendChild(addCardWithFilters('Unique Users Per Day', 'UniqueUsersPerDayOnDateAndUsersAndEvents', (t: DG.DataFrame) => DG.Viewer.lineChart(t).root));
            errorListCard.appendChild(addCardUsingDataframe('Error Types List', errors, (t: DG.DataFrame) => DG.Viewer.grid(t).root));
            
            usageCard.appendChild(
                    addCardWithFilters('Usage', 'EventsOnDateAndUsersAndEvents', (t: DG.DataFrame) => {
                    subscribeOnTableWithEvents(t);
                    return DG.Viewer.scatterPlot(t, {'color': 'user'}).root;
                    })
            );

            eventTypeCard.appendChild(addCardWithFilters('Event Types', 'EventsSummaryOnDate', (t: DG.DataFrame) => DG.Viewer.barChart(t, {valueAggrType: 'avg'}).root));
            errorTypeCard.appendChild(addCardUsingDataframe('Error Types', errors, (t: DG.DataFrame) => DG.Viewer.barChart(t, {valueAggrType: 'avg'}).root));
            testTrackingCard.appendChild(addCardWithFilters('Test Tracking', 'ManualActivityByDate', (t: DG.DataFrame) => DG.Viewer.grid(t).root));
            
            this.root.append(ui.divV([
                ui.divH([usersSummaryCard,eventsSummaryCard,errorsSummaryCard]),
                ui.divH([uniqueUsersCard,uniqueUsersPerDayCard]),
                ui.divH([usageCard]),
                ui.divH([eventTypeCard,errorTypeCard]),
                ui.divH([errorListCard]),
                ui.divH([testTrackingCard])
            ])
            );
            ui.setUpdateIndicator(this.root, false);
        }

        (async () => {
            await showUsage();
        })();

        this.toolbox = accToolbox.root;
        
        // return html card
        function addCardUsingDataframe(cardName:string, dataframe:DG.DataFrame, f:any, supportUsers = true) {
            let host = ui.block([],'d4-item-card card');
            host.appendChild(ui.h1(cardName));

            if (cardName === 'Errors')
            grok.data.detectSemanticTypes(dataframe);
            host.appendChild(f(dataframe));
            return host;
        }

        // return current filter
        function getCurrentFilter() {
            return {'date': dateInput.value, 'events': [eventsInput.value], 'isExactly': isExactlyInput.value, 'users': [usersInput.value]};
        }
        
        // return html card with aplied filters
        function addCardWithFilters(cardName:string, queryName:string, viewer:any) {
            let host = ui.block([],'d4-item-card card');
            host.appendChild(ui.h1(cardName));
            let loader = ui.loader();
            host.appendChild(loader);
    
            let filter = getCurrentFilter();
    
            grok.data.query('UsageAnalysis:' + queryName, filter).then((t) => {
              if (cardName === 'Errors')
                grok.data.detectSemanticTypes(t);
              host.appendChild(viewer(t));
              host.removeChild(loader);
            });
            return host;
          }

        //  
        function getUniqueUsers(cardName:string){
            let host = ui.block([],'d4-item-card card');
            host.style.overflow="auto";
            host.style.height="94.5%";
            host.appendChild(ui.h1(cardName));
            let loader = ui.loader();
            host.appendChild(loader);
            grok.data.query('UsageAnalysis:UniqueUsersByDate', {'date': dateInput.value})
            .then((t: DG.DataFrame) => {
                let ids = Array.from(t.getCol('id').values());
                grok.dapi.getEntities(ids).then((users) => {
                host.appendChild(ui.list(users));
                });
            });
            host.removeChild(loader);
            return host;
        }

        
    }

    debounce(fn:any, time:number) {
        let timeout:any;
        return  () => {
          const functionCall = () => fn.apply(this, arguments);
          clearTimeout(timeout);
          timeout = setTimeout(functionCall, time);
        }
    }

}

