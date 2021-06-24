# """Testing type hinting"""

from typing import TypeVar

T = TypeVar('T', bound='Foo')

class Foo:
    def __init__(self, value: int) -> None:
        self.value = value

    def new_with_value(self: T, value: int) -> T:
        return self.__class__(value)

    def increment_value(self: T) -> T:
        return self.new_with_value(self.value + 1)

class Bar(Foo):
    def __init__(self, value: int) -> None:
        self.value = value


def return_a_new_one(x: T) -> T:
    return x.new_with_value(x.value + 1)

f = Bar(1)
g = return_a_new_one(f)
print(g.value)
print(type(g))
# from typing import TypeVar

# T = TypeVar('T', bound='Shape')

# class Shape:
#     def set_scale(self: T, scale: float) -> T:
#         self.scale = scale
#         return self

# class Circle(Shape):
#     def set_radius(self, r: float) -> 'Circle':
#         self.radius = r
#         return self

# class Square(Shape):
#     def set_width(self, w: float) -> 'Square':
#         self.width = w
#         return self

# circle = Circle().set_scale(0.5).set_radius(2.7)  # type: Circle
# square = Square().set_scale(0.5).set_width(3.2)  # type: Square


x = [1,2,3,4,5,6]

for i, val1 in enumerate(x):
    for val2 in x[i:]:
        print(val1, val2)



