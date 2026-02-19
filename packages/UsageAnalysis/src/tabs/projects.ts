import {UaView} from "./ua";
import {UaToolbox} from "../ua-toolbox";
import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import {UaFilterableQueryViewer} from "../viewers/ua-filterable-query-viewer";
import {loadUsers, setupUserIconRenderer} from "../utils";

export class ProjectsView extends UaView {
    constructor(uaToolbox?: UaToolbox) {
        super(uaToolbox);
        this.name = 'Projects';
    }

    async initViewers(path?: string): Promise<void> {
        const users = await loadUsers();
        let syncing = false;

        const syncViewers = (projectName: string, source: 'scatter' | 'bar' | 'grid') => {
            if (syncing) return;
            syncing = true;
            try {
                if (source !== 'scatter') {
                    const df = projectsViewer.viewer?.dataFrame;
                    if (df) {
                        const col = df.getCol('project_name');
                        df.selection.init((i) => col.get(i) === projectName);
                    }
                }
                if (source !== 'bar') {
                    const df = uniqueUsersPerProject.viewer?.dataFrame;
                    if (df) {
                        const col = df.getCol('project_name');
                        df.selection.init((i) => col.get(i) === projectName);
                    }
                }
                if (source !== 'grid') {
                    const df = dailyProjectAccess.viewer?.dataFrame;
                    if (df) {
                        const col = df.getCol('project_name');
                        df.selection.init((i) => col.get(i) === projectName);
                    }
                }
            }
            finally {
                syncing = false;
            }
        };

        const projectsViewer = new UaFilterableQueryViewer({
            filterSubscription: this.uaToolbox.filterStream,
            name: 'Projects Usage',
            queryName: 'AccessCountPerPeriodPerProject',
            processDataFrame: (t: DG.DataFrame) => {
                t.onCurrentRowChanged.subscribe(() => {
                    if (t.currentRowIdx < 0) return;
                    syncViewers(t.get('project_name', t.currentRowIdx), 'scatter');
                });
                return t;
            },
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
                const viewer = DG.Viewer.barChart(t, {
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
                viewer.onEvent('d4-bar-chart-on-category-clicked').subscribe((args: any) => {
                    syncViewers(args.args.options.categories[0], 'bar');
                });
                return viewer;
            }
        });

        const dailyProjectAccess = new UaFilterableQueryViewer({
            filterSubscription: this.uaToolbox.filterStream,
            name: 'Daily Project Access',
            queryName: 'DailyProjectAccess',
            processDataFrame: (t: DG.DataFrame) => {
                t.getCol('project_name').setTag('friendlyName', 'Project Name');
                t.getCol('access_date').setTag('friendlyName', 'Date');
                t.getCol('user_name').setTag('friendlyName', 'User');
                t.onCurrentRowChanged.subscribe(() => {
                    if (t.currentRowIdx < 0) return;
                    syncViewers(t.get('project_name', t.currentRowIdx), 'grid');
                });
                return t;
            },
            createViewer: (t: DG.DataFrame) => {
                const grid = DG.Viewer.grid(t, {
                    'title': 'Daily Project Access',
                    'description': 'Shows which unique users opened each project per day.',
                    'descriptionVisibilityMode': DG.VisibilityMode.Never
                });
                grid.sort(['access_date'], [false]);
                setupUserIconRenderer(grid, users, ['user_name']);
                return grid;
            }
        });

        this.viewers.push(projectsViewer);
        this.viewers.push(uniqueUsersPerProject);
        this.viewers.push(dailyProjectAccess);
        this.root.append(projectsViewer.root);

        this.root.append(
            ui.splitV([
                    projectsViewer.root,
                    ui.splitH([
                        ui.box(uniqueUsersPerProject.root, {style: {height: '100%'}}),
                        ui.box(dailyProjectAccess.root, {style: {height: '100%'}}),
                    ])
                ]
            )
        );
        // this.root.append(packagesTimeViewer.root);
    }
}