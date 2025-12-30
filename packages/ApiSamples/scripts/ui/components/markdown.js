// Render markdown 

let document = `# Intro
Datagrok unlocks the value of the complex data by empowering non-technical users to: \n 
  * [Discover](/help/discover/fair). 
  * [Cleanse](/help/transform/data-wrangling). 
  * [Visualize](/help/visualize/viewers).
  * [Explore data](/help/#explore).
  * [Build and deploy predictive models](/help/learn/predictive-modeling).
  * [Collaborate with others](/help/collaborate/sharing).
`;

grok.shell.newView('Markdown example', [
  ui.markdown(document)
]);