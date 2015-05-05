from multiprocessing import Process, Value, Lock

class Counter():
    """Thread safe counter"""

    def __init__(self, initval=0, typ='i'):
        """Constructor"""
        self.val = Value(typ, initval)
        self.lock = Lock()

    def increment(self):
        """Increments counter"""
        with self.lock:
            self.val.value += 1

    def value(self):
        """Returns value"""
        with self.lock:
            return self.val.value

    def add(self, val):
        """Adds to counter"""
        with self.lock:
            self.val.value += val

    def __str__(self):
        """String representation"""
        return str(self.val.value)
