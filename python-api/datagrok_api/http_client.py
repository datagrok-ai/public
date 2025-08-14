import requests
from typing import Any, Dict, Optional

class HttpClient:
    def __init__(self, base_url: str, api_key: str):
        self.session = requests.Session()
        self.base_url = base_url
        self.session.headers.update({
            'Authorization': api_key
        })

    def request(self, method: str, endpoint: str, 
                headers: Optional[Dict[str, str]] = None, **kwargs: Any) -> requests.Response:
        """
        Make an HTTP request to the Datagrok API.

        Parameters
        ----------
        method : str
            HTTP method to use (e.g., "GET", "POST", "PUT", "DELETE").
        endpoint : str
            API endpoint path, relative to the base URL.
        headers : dict of str to str, optional
            Optional headers to include or override for this request.
        **kwargs : Any
            Additional arguments to pass to `requests.request`
            (e.g., `params`, `json`, `data`, `files`, etc.).

        Returns
        -------
        requests.Response
            The response returned by the API.

        Raises
        ------
        Exception
            If the API returns an error or if the response is not successful
        """

        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        merged_headers = {
            'Content-Type': 'application/json',
            **self.session.headers,
            **(headers or {}),
        }
        response = self.session.request(method, url, headers=merged_headers, **kwargs)
        response.raise_for_status()
        if 'api-error' in response.headers:
            content = response.content.decode()
            if 'ApiError' in content:
                raise Exception(content)
            else:
                raise Exception("Something went wrong during request") 
        return response
    
    def get(self, endpoint: str, 
                headers: Optional[Dict[str, str]] = None, **kwargs: Any):
        return self.request('GET', endpoint, headers=headers, **kwargs)
    
    def post(self, endpoint: str, 
                headers: Optional[Dict[str, str]] = None, **kwargs: Any):
        return self.request('POST', endpoint, headers=headers, **kwargs)
    
    def delete(self, endpoint: str, 
                headers: Optional[Dict[str, str]] = None, **kwargs: Any):
        return self.request('DELETE', endpoint, headers=headers, **kwargs)

    def close(self):
        self.session.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
    