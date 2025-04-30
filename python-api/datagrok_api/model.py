import uuid

class Model:
    def __init__(self, id):
        self.id = id

    def ensure_id(self) -> str:
        if self.id is None:
            self.id = str(uuid.uuid1())
        return self.id    