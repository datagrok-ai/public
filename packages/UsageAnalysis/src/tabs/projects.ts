import {UaView} from "./ua";
import {UaToolbox} from "../ua-toolbox";
import * as ui from "datagrok-api/ui";
import * as grok from "datagrok-api/grok";
import * as DG from "datagrok-api/dg";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";

export class ProjectsView extends UaView {
    constructor(uaToolbox: UaToolbox) {
        super(uaToolbox);
        this.name = 'Projects';
    }

    async initViewers(path?: string): Promise<void> {
        const projectsViewer = new UaFilterableQueryViewer({
            filterSubscription: this.uaToolbox.filterStream,
            name: 'Projects Usage',
            queryName: 'AccessCountPerPeriodPerProject',
            createViewer: (t: DG.DataFrame) => {
                return DG.Viewer.scatterPlot(t, {
                    x: 'period',
                    y: 'project_name',
                    size: 'access_count',
                    color: 'project_name',
                    jitterSize: 5,
                    markerMinSize: 10,
                    markerMaxSize: 30,
                    showColorSelector: false,
                    showSizeSelector: false,
                    showXSelector: false,
                    showYSelector: false,
                    invertYAxis: true,
                });
            }});

        const uniqueUsersPerProject = new UaFilterableQueryViewer({
            filterSubscription: this.uaToolbox.filterStream,
            name: 'Unique Users Per Project',
            queryName: 'UniqueUsersPerProject',
            createViewer: (t: DG.DataFrame) => {
                return DG.Viewer.barChart(t, {
                    'valueColumnName': 'unique_users',
                    'valueAggrType': 'sum',
                    'barSortType': 'by value',
                    'barSortOrder': 'desc',
                    'showValueAxis': false,
                    'showValueSelector': false,
                    'splitColumnName': 'project_name',
                    'showCategoryValues': true,
                    'showCategorySelector': false,
                    'stackColumnName': '',
                    'showStackSelector': false,
                    'title': 'Unique Users Per Project',
                    'rowSource': DG.RowSet.All,
                    'description': 'Counts the number of distinct users accessing each project.',
                    'descriptionVisibilityMode': DG.VisibilityMode.Never
                });
            }
        });

        const accessFrequencyPerProject = new UaFilterableQueryViewer({
            filterSubscription: this.uaToolbox.filterStream,
            name: 'Access Frequency',
            queryName: 'UserAccessFrequencyPerProject',
            processDataFrame: (t: DG.DataFrame) => {
                t.getCol('project_name').setTag('friendlyName', 'Project Name');
                t.getCol('unique_users').setTag('friendlyName', 'Avg Unique Users Daily');
                t.getCol('avg_days_between_accesses').setTag('friendlyName', 'Days Between Access');
                return t;
            },
            createViewer: (t: DG.DataFrame) => {
                const grid = DG.Viewer.grid(t, {
                    'title': 'Access Frequency Daily',
                    'description': 'Provides average number of days between user accesses and average users count daily for projects that have more than one access.',
                    'descriptionVisibilityMode': DG.VisibilityMode.Never
                });
                grid.sort(['unique_users', 'avg_days_between_accesses'], [false, true]);
                return grid;
            }
        });

        this.viewers.push(projectsViewer);
        this.viewers.push(uniqueUsersPerProject);
        this.viewers.push(accessFrequencyPerProject);
        // this.viewers.push(packagesTimeViewer);
        // packagesTimeViewer.root.style.display = 'none';
        this.root.append(projectsViewer.root);

        this.root.append(
            ui.splitV([
                    projectsViewer.root,
                    ui.splitH([
                        ui.box(uniqueUsersPerProject.root, {style: {height: '100%'}}),
                        ui.box(accessFrequencyPerProject.root, {style: {height: '100%'}}),
                    ])
                ]
            )
        );
        // this.root.append(packagesTimeViewer.root);
    }
}