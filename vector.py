from math import sqrt, acos, pi, radians
from decimal import Decimal, getcontext

getcontext().prec = 30


class Vector(object):
	CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector.'
	def __init__(self, coordinates):
		try:
			if not coordinates:
				raise ValueError
			self.coordinates = tuple([Decimal(x) for x in coordinates])
			self.dimension = len(coordinates)
		except ValueError:
			raise ValueError("The coordinates must be nonempty.")
		except TypeError:
			raise TypeError("The coordinates must be iterable.")
	def __str__(self):
		"""Returns input as a string"""
		return "Vector: {}".format(self.coordinates)
	def __eq__(self, v):
		"""This method allows you to compare two instances for equality."""
		return self.coordinates == v.coordinates
	def plus(self, v):
		return [x+y for x,y in zip(self.coordinates, v.coordinates)]
	def minus(self, v):
		return [x-y for x,y in zip(self.coordinates, v.coordinates)]
	def scalar_multiply(self, c):
		scaled_coordinates = [c*x for x in self.coordinates]
		return Vector(scaled_coordinates)
	def magnitude(self):
		return sqrt(sum([x**2 for x in self.coordinates]))
	def normalized(self):
		try:
			magnitude = Decimal(self.magnitude())
			return self.scalar_multiply(Decimal(1.0)/magnitude)
		except ZeroDivisionError:
			raise Exception("The vector cannot be the zero vector.")
	def dot_product(self, v):
		return sum((x*y for x,y in zip(self.coordinates, v.coordinates)))
	def angle(self, v, in_degrees=False):
		try:
			v1 = self.normalized()
			v2 = v.normalized()
			degrees = acos(round(v1.dot_product(v2), 12))
			if in_degrees:
				return radians(degrees)
			else:
				return degrees
		except Exception as e:
			if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
				raise Exception('Cannot comput an angle with the zero vector.')
			else:
				raise e
	def is_zero(self, tolerance=1e-10):
		return self.magnitude() < tolerance
	def is_parallel_to(self, v):
		return (self.is_zero() or
				v.is_zero() or
				self.angle(v) == 0 or
				self.angle(v) == pi)
	def is_orthogonal_to(self, v, tolerance=1e-10):
		return abs(self.dot_product(v)) < tolerance
	def component_orthogonal_to(self, basis):
		"""Find v perp"""
		try:
			v_parallel = self.component_parallel_to(basis)
			return self.minus(v_parallel)
		except Exception as e:
			if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
				raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
			else:
				raise e
	def component_parallel_to(self, basis):
		"""Find v parallel"""
		try:
			b = basis.normalized()
			scalar = v1.dot_product(b)
			return b.scalar_multiply(scalar)
		except Exception as e:
			if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
				raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
			else:
				raise e
	def cross_product(self, v):
		"""Assuming three dimensions"""
		x_1, y_1, z_1 = self.coordinates
		x_2, y_2, z_2 = v.coordinates
		new_coordinates = [	y_1*z_2 - y_2*z_1, -(x_1*z_2 - x_2*z_1), x_1*y_2 - x_2*y_1]
		return Vector(new_coordinates)
	def area_parallelogram(self, base):
		"""Assuming three dimensions"""
		height_vector = self.cross_product(base)
		area = height_vector.magnitude()
		return area
	def area_triangle(self, base):
		"""Assuming three dimensions"""
		return self.area_parallelogram(base) * 0.5