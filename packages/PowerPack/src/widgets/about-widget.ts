/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export class AboutWidget extends DG.Widget {
  caption: string;

  constructor() {
    super(ui.box());

    // properties
    this.caption = super.addProperty('caption', DG.TYPE.STRING, 'Why Datagrok?');

    let tabControl = ui.tabControl({
    'PLATFORM' : () => ui.panel([ui.markdown(intro)]),
    'USERS' : () => ui.panel([ui.markdown(introUsers)]),
    'IT' : () => ui.panel([ui.markdown(introIT)]),
    'DATA SCIENCE' : () => ui.panel([ui.markdown(introDS)]),
    'CHEMISTRY' : () => ui.panel([
      ui.markdown(introChem),
      ui.div([
        ui.link('Go ahead, sketch a molecule!', async () => 
        // @ts-ignore 
        ui.wait(()=>{
          grok.functions.call('CmdSketcher')
          //cmdChemSketcher.run()
        })
      )], {style:{marginTop:'10px'}})
    ]),
    'DEVELOPMENT' : () => ui.panel([ui.markdown(introDev)]),
    'BUSINESS' : () => ui.panel([ui.markdown(introLeaders)]),
    //
    })
    let currentPane = 0;
    let tabControlButtonLeft = ui.div([
      ui.button(ui.iconFA('chevron-left'), ()=>{
        if (currentPane != 0){
          currentPane--;
          tabControl.currentPane = tabControl.panes[currentPane];
        }
      })]);
      let tabControlButtonRight = ui.div([
      ui.button(ui.iconFA('chevron-right'), ()=>{
        if (currentPane<6){
          currentPane++;
          tabControl.currentPane = tabControl.panes[currentPane];
        }
      })
    ]);
    tabControlButtonRight.style.float = 'Right';
    tabControlButtonRight.style.margin = '0px 10px 0px 0px';
    tabControlButtonRight.style.color = 'var(--blue-1)';
    tabControlButtonLeft.style.float = 'Left';
    tabControlButtonLeft.style.margin = '0px 0px 0px 10px';
    tabControlButtonLeft.style.color = 'var(--blue-1)';

    tabControl.header.className = 'd4-tab-header-stripe about-widget';
    tabControl.header.prepend(tabControlButtonLeft);
    tabControl.header.appendChild(tabControlButtonRight);
    //tabControl.root.className = 'd4-tab-host d4-tab-vertical';
    this.root.appendChild(tabControl.root);
  }
}

const intro =  `Datagrok unlocks the value of the complex data by empowering 
non-technical users to 
[discover](/help/discover/fair), 
[cleanse](/help/transform/data-wrangling), 
[visualize](/help/visualize/viewers), 
[explore](/help/#explore) data, 
build and deploy [predictive models](/help/learn/predictive-modeling), and collaborate with others. \n
Datagrok's [performance](/help/develop/performance) is unmatched 
due to our technology that allows to work with tens of millions of data points interactively
right in the web browser. The platform is flexible and extensible, and integrates easily with existing systems. If anything
is missing, extend the platform with  [plugins and applications](/help/develop/develop). \n
The platform is a lot more than a tool for data analysis - it is an ecosystem designed
for users, IT, data scientists, chemists, and developers not only to collaborate, but to 
create value for other types of users as well. Click on the tabs to the left to see how 
the platform would benefit different users.`;

const introUsers = `Access any data source, anytime. Use the same platform for
data discovery, self-service data analytics, creating and sharing dashboards.
Run fit-for-purpose applications developed on top of the platform.
* Keep your files in order, access from any computer, and share with others
* Get [context-specific suggestions](/help/discover/data-augmentation) for analyzing data in front of you`;

const introIT = `Software is eating the world, and IT is a critical competence for all enterprises.
Datagrok enables IT to become a partner with business by solving business problems, at the same time 
providing the following:
* Enterprise-grade [security](/help/govern/security)
* Central management of database connections and queries
* Integration with the existing systems 
* Data governance, data lineage, and [audit](/help/govern/audit)
* Lightweight, managed deployment of scientific methods and predictive models
* [Usage analysis](/help/govern/usage-analysis)`;

const introDS = `
Datagrok lets you easily connect to data, visualize and explore it, build predictive models, 
and deploy your work to production at a speed previously unheard of. Unlike other platforms,
"deployment" in Datagrok is not simply exposing an endpoint API. \n 
Instead, it is an integrated 
experience, involving data discovery, visualization, and augmentation. See our
[YouTube channel](https://www.youtube.com/watch?v=tVwpRB8fikQ) where we would start from scratch
and end up training a predictive model and producing a rich application using it.
* Interactively retrieve, clean, and transform datasets
* Enrich db query results with the augmented data, and give them to the end users 
* Unmatched capabilities for [exploratory data analysis](/help/explore/exploratory-data-analysis)
* Deploy [R, Python, Julia](/help/compute/scripting) scripts into production''';`

const introChem = `''Datagrok provides [first-class support for small molecules](/help/domains/chem/cheminformatics).
* [Chemically aware viewers](/help/domains/chem/chemically-aware-viewers)
* Rich set of [molecular descriptors](/help/domains/chem/descriptors)
* Support for RDKit or OpenChemLib  
* [Sketch molecules](/help/domains/chem/sketcher) 
* [R-group analysis](/help/domains/chem/r-group-analysis)
* Chemical [Similarity](/help/domains/chem/similarity-search) and [diversity](/help/domains/chem/diversity-search) analyses`;

const introDev = `
Think about Datagrok as an OS for data. Just like the way Linux enables application developers to be productive
by handling low-level stuff and introducing high-level concepts like processes, pipes, and files, Datagrok
does the same for data. By leveraging our proprietary high-performance [data engine](/help/develop/performance), we provide
efficient solutions for [data connectivity](/help/access/databases#database-manager), [visualization](/help/visualize/viewers), 
[scientific computations](/help/compute/scripting), and
[predictive modeling](/help/learn/predictive-modeling). \n
Extend the platform by [writing functions](/help/compute/scripting) in SQL, R, Python, or Julia. 
This lightweight process will take you a long way; some applications that take months to develop could 
now be [done in a matter of days](https://www.youtube.com/watch?v=tVwpRB8fikQ). \n
For deeper integration and customization, use our [JS API](/help/develop/develop) to take full control
of the platform.`;

const introLeaders = `Connect business areas with IT, data science, and application developers.
* Take control of your data and analytics landscape by connecting with existing systems
* Break organizational barriers by combining data and algorithms from disparate data silos and functional areas
* Start small, and organically grow your ecosystem
* Keep your hand on the pulse, and get insights from usage analysis
* Leverage cross-domain competencies and capabilities`;
