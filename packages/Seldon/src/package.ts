/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();

//name: SeldonApply
//top-menu: Seldon|Apply
export async function SeldonApply() {

  // 1. Obtain the list of deployments
  let p = new DG.Package();
  p.name = 'Seldon';
  const credentialsResponse = await p.getCredentials();
  const parameters : {[index: string] : any} = credentialsResponse.parameters;
  const seldonUser : string = parameters['seldonUser'];
  const seldonPassword : string = parameters['seldonPassword'];
  const seldonHost : string = parameters['seldonHost'];
  const seldonOIDCServer : string = parameters['seldonOIDCServer'];
  const seldonClientID : string = parameters['seldonClientID'];
  const seldonNamespace : string = parameters['seldonNamespace']; // get_cluster
  let seldonDeploymentsDf = await grok.functions.call(
    "Seldon:SeldonGetDeploymentsByNamespacePy", {
      'seldonUser': seldonUser,
      'seldonPassword': seldonPassword,
      'seldonHost': seldonHost,
      'seldonOIDCServer': seldonOIDCServer,
      'seldonClientID': seldonClientID,
      'seldonNamespace': seldonNamespace
  });

  const seldonDeploymentsList = seldonDeploymentsDf.columns.byName('seldonDeployments').toList();
  let choiceInputDeployment = ui.choiceInput('Model', seldonDeploymentsList[0], seldonDeploymentsList);
  ui.dialog({
    title: 'Seldon models'
  })
  .add(ui.span([
    `Select a model from the namespace '${seldonNamespace}'`
  ]))
  .add(ui.div([choiceInputDeployment]))
  .onOK(async () => {
    let seldonInputDf = grok.shell.tables[0]; // sic!
    // 2. Apply the selected deployment model
    
    
    let f = await grok.functions.eval("Seldon:SeldonApplyModelPy");
    let call = f.prepare({      
        'seldonUser': seldonUser,
        'seldonPassword': seldonPassword,
        'seldonHost': seldonHost,
        'seldonOIDCServer': seldonOIDCServer,
        'seldonClientID': seldonClientID,
        'seldonNamespace': seldonNamespace,
        'seldonDeploymentName': choiceInputDeployment.value,
        'seldonInput': seldonInputDf
    });
    await call.call();
    const seldonErrorState = call.getParamValue('seldonErrorState');
    if (seldonErrorState === '') {
      const seldonResultDf = call.getParamValue('seldonResult');
      for (let col of seldonResultDf.columns)
        seldonInputDf.columns.insert(col);
    }
  })
  .show();

}
