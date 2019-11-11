class GrokApiExamplesPackage extends GrokPackage {

    // Guaranteed to get called exactly once before the execution of any function below
    init() { console.log('Examples package initialized.'); }
    
    //description: test
    test() {
       alert('Examples package test');
    }

}
