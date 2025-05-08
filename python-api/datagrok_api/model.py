import uuid

class Model:
    """Base class for all Datagrok models.
    
    This class provides basic functionality for all model classes in the Datagrok API.
    It handles unique identification of model instances and ensures each instance
    has a valid ID.
    
    Attributes
    ----------
    id : str
        Unique identifier for the model instance. Can be None for new instances
        that haven't been saved to the server yet.
    """
    
    def __init__(self, id):
        """Initialize a new Model instance.
        
        Parameters
        ----------
        id : str
            Unique identifier for the model instance. Can be None for new instances.
        """
        self.id = id

    def ensure_id(self) -> str:
        """Ensure the model has a valid ID, generating one if necessary.
        
        If the model's ID is None, generates a new UUID v1 and assigns it to the model.
        This is typically used before saving a new model instance to the server.
        
        Returns
        -------
        str
            The model's ID (either existing or newly generated)
        """
        if self.id is None:
            self.id = str(uuid.uuid1())
        return self.id    