from multiprocessing import Process, Value, Lock

class Counter():
    def __init__(self, initval=0, typ='i'):
        self.val = Value(typ, initval)
        self.lock = Lock()

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value

    def add(self, val):
        with self.lock:
            self.val.value += val

    def __str__(self):
        return str(self.val.value)
