"""."""

from typing import List


# A custom type hint about "Vectors", which are just lists containing floats.
Vector = List[float]


def vector_add(v: Vector, w: Vector) -> Vector:
    """Adds elements of a vector to produce a new vector."""
    assert len(v) == len(w), "Vectors must be the same length"
    return [v_i + w_i for (v_i, w_i) in zip(v, w)]


assert vector_add([1, 2, 3], [4, 5, 6]) == [5, 7, 9]


def vector_subtract(v: Vector, w: Vector) -> Vector:
    """Subtracts elements of a vector to produce a new vector."""
    assert len(v) == len(w), "Vectors must be the same length"
    return [v_i - w_i for (v_i, w_i) in zip(v, w)]


assert vector_subtract([4, 5, 6], [1, 2, 3]) == [3, 3, 3]
assert vector_subtract([5, 7, 9], [4, 5, 6]) == [1, 2, 3]


def vector_sum(vectors: List[Vector]) -> Vector:
    """Sums elements of all input vectors."""
    assert vectors, "No vectors provided!"
    num_elements = len(vectors[0])
    assert all(len(v) == num_elements for v in vectors), "Vectors are different sizes!"
    return [sum(vector[i] for vector in vectors) for i in range(num_elements)]


assert vector_sum([[1, 2], [3, 4], [5, 6], [7, 8]]) == [16, 20]


def scalar_multiply(c: float, v: Vector) -> Vector:
    """Multiply a vector by a scalar."""
    return [c * v_i for v_i in v]


assert scalar_multiply(3, [1, 1, 1]) == [3, 3, 3]


def vector_mean(vectors: List[Vector]) -> Vector:
    """Computes the element-wise mean of a set of vectors."""
    n = len(vectors)
    return scalar_multiply(1 / n, vector_sum(vectors))


assert vector_mean([[1, 2], [3, 4], [5, 6]]) == [3, 4]


# Sum of componentwise products.
def dot(v: Vector, w: Vector) -> float:
    """Computes the dot product of two vectors."""
    assert len(v) == len(w), "Vectors must be the same length!"
    return sum(v_i * w_i for v_i, w_i in zip(v, w))


assert dot([1, 2, 3], [4, 5, 6]) == 32


def sum_of_squares(v: Vector) -> float:
    """Computes the sum of every element squared."""
    return dot(v, v)


assert sum_of_squares([1, 2, 3]) == 14
