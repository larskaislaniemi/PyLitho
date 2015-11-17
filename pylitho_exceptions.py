class PyLerr_OutsideRange(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PyLerr_Undefined(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PyLerr_TypeError:
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
