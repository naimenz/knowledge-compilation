def f(*args: int):
    for arg in args:
        print(arg)
    return args

print(f(1,2,3))

