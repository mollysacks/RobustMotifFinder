import time

def ReportRuntime(f):
    def g(*args):
        t0 = time.time()
        output = f(*args)
        t1 = time.time()
        if output:
            return *output, t1-t0
        else:
            return t1-t0
    g.__name__ = f.__name__
    return g

@ReportRuntime
def Test(a, b, c):
    print(f"inputs are {a} {b} {c}")
    return a, b, c


Test(1, 2, 3)
