class ShareResponse:
    """
    Represents the result of a share operation in the Datagrok API.

    This class encapsulates the status of a share request, including 
    successful shares, already shared items, and failures.

    Parameters
    ----------
    data : dict
        Dictionary containing share response data from the Datagrok API.
        Expected keys:
            - 'status' (str): The status of the operation, one of:
                * 'success' — All items were shared successfully.
                * 'partial_success' — Some items were shared successfully, some failed.
                * 'failed' — No items were shared.
            - 'shared' (list): Items that were successfully shared.
            - 'alreadyShared' (list): Items that were already shared before the request.
            - 'failed' (list): Items that failed to be shared.

    Attributes
    ----------
    status : str
        The overall status of the share operation.
    shared : list
        Items successfully shared.
    already_shared : list
        Items already shared before the request.
    failed : list
        Items that failed to be shared.

    Methods
    -------
    is_success() -> bool
        Returns True if all items were successfully shared.
    is_partial() -> bool
        Returns True if the operation had both successes and failures.
    is_failed() -> bool
        Returns True if the operation failed entirely.
    raise_for_failure() -> None
        Raises a RuntimeError if the operation failed.
    """
    def __init__(self, data: dict):
        self.status = data.get('status')
        self.shared = data.get('shared', [])
        self.already_shared = data.get('alreadyShared', [])
        self.failed = data.get('failed', [])

    def is_success(self):
        return self.status == 'success'

    def is_partial(self):
        return self.status == 'partial_success'

    def is_failed(self):
        return self.status == 'failed'

    def raise_for_failure(self):
        if self.is_failed():
            raise RuntimeError(f"Sharing failed: {self.failed}")
        