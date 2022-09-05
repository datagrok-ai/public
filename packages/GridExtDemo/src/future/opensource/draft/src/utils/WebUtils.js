export class WebUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}


WebUtils.makeGetRequest = async function(url) {
    const myHeaders = new Headers();
    myHeaders.append('Content-Type', 'application/json');
    const request = new Request(url, {method: 'GET', credentials : "include", headers: myHeaders});

    let response = null;
    try{response = await fetch(request);}
    catch (error)
    {
        console.log(error.message);
        throw error;
    }


    if (response.ok)
    {
        const obJSON = await response.json();
        const strJSon = JSON.stringify(obJSON);
        console.log(strJSon);

        return obJSON;
    }
    else
    {
        let bRedirected = response.redirected;
        let urlResp =  response.url;

        console.log("Something went wrong: Response Code: " + response.status + " " + response.error);
        return null;
    }
}

WebUtils.makeWebRequest = async function(url, body) {

    const myHeaders = new Headers();
    myHeaders.append('Content-Type', 'application/json');
    const request = new Request(url, {method: 'POST', credentials : "include", headers: myHeaders, body: body});

    let response = null;
    try{response = await fetch(request);}
    catch (error)
    {
        console.log(error.message);
        throw error;
    }


   if (response.ok)
   {
    const obJSON = await response.json();
    const strJSon = JSON.stringify(obJSON);
    //console.log(strJSon);

    return obJSON;
   }
   else
   {
    let bRedirected = response.redirected;
    let urlResp =  response.url;

    console.log("Something went wrong: Response Code: " + response.status + " " + response.error);
    return null;
   }
}
