from typing import Optional
import uuid
from abc import ABC, abstractmethod
from datetime import datetime

class Model(ABC):
    """Base class for all Datagrok models.
    
    This class provides basic functionality for all model classes in the Datagrok API.
    It handles unique identification of model instances and ensures each instance
    has a valid ID.
    
    Parameters
    ----------
    id : str
        Unique identifier for the model instance. Can be None for new instances
        that haven't been saved to the server yet.
    """
    
    def __init__(self, id: Optional[str] = None):
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


class NamedModel(Model, ABC):
    """NamedModel class that contains some common fields of Datagrok entities.

    Parameters
    ----------
    name : Optional[str]
        Internal name of the model, used as an identifier.
    id : Optional[str]
        Unique identifier for the model instance. Can be None for new instances.
    friendly_name : Optional[str]
        Human-readable name for display purposes. Defaults to None.
    created_on : Optional[datetime]
        Timestamp indicating when the model was created. Defaults to None.
    updated_on : Optional[datetime]
        Timestamp indicating when the model was last updated. Defaults to None.
    namespace : Optional[str]
        Optional namespace prefix to distinguish models with the same name. Defaults to None.
    """
    def __init__(self, name: Optional[str] = None,
                 id: Optional[str] = None, 
                 friendly_name: Optional[str] = None,
                 created_on: Optional[datetime] = None,
                 updated_on: Optional[datetime] = None,
                 namespace: Optional[str] = None):
        super().__init__(id)
        self.name = name
        self.friendly_name = friendly_name
        self.created_on = created_on
        self.updated_on = updated_on
        self.namespace = namespace

    @property
    def grok_name(self):
        return (self.namespace or '') + (self.name or '')  
    
    @abstractmethod
    def to_dict(self) -> dict:
        """Convert the model to a dictionary representation."""
        pass    
