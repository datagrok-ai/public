class ShareResponse:
    """
    Represents the result of a share operation in the Datagrok API.

    This class encapsulates the status of a share request, including 
    successful shares, already shared items, and failures.

    Attributes
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
    status : str
        The overall status of the share operation.
    shared : list
        Items successfully shared.
    already_shared : list
        Items already shared before the request.
    failed : list
        Items that failed to be shared.
    """
    def __init__(self, data: dict):
        self.status = data.get('status')
        self.shared = data.get('shared', [])
        self.already_shared = data.get('alreadyShared', [])
        self.failed = data.get('failed', [])

    def is_success(self):
        """
        Returns True if all items were successfully shared.
        """
        return self.status == 'success'

    def is_partial(self):
        """
        Returns True if the operation had both successes and failures.
        """
        return self.status == 'partial_success'

    def is_failed(self):
        """
         Returns True if the operation failed entirely
        """
        return self.status == 'failed'

    def raise_for_failure(self):
        """
        Raises a RuntimeError if the operation failed.
        """
        if self.is_failed():
            raise RuntimeError(f"Sharing failed: {self.failed}")
        