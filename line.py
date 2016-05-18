from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30

class Line(object):
	NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
	def __init__(self, normal_vector=None, constant_term=None):
		self.dimension = 2
		if not normal_vector:
			all_zeros = ['0']*self.dimension
			normal_vector = Vector(all_zeros)
		self.normal_vector = normal_vector.coordinates
		if not constant_term:
			constant_term = Decimal('0')
		self.constant_term = Decimal(constant_term)
		self.set_basepoint()
	def set_basepoint(self):
		try:
			n = self.normal_vector
			c = self.constant_term
			basepoint_coords = ['0']*self.dimension
			initial_index = Line.first_nonzero_index(n)
			initial_coefficient = n[initial_index]
			basepoint_coords[initial_index] = c / initial_coefficient
			self.basepoint = Vector(basepoint_c2oords)
		except Exception as e:
			if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
				self.basepoint = None
			else:
				raise e
	def __str__(self):
		num_decimal_places = 3
		def write_coefficient(coefficient, is_initial_term=False):
			coefficient = round(coefficient, num_decimal_places)
			if coefficient % 1 == 0:
				coefficient = int(coefficient)
			output = ''
			if coefficient < 0:
				output += '-'
			if coefficient > 0 and not is_initial_term:
				output += '+'
			if not is_initial_term:
				output += ' '
			if abs(coefficient) != 1:
				output += '{}'.format(abs(coefficient))
			return output
		n = self.normal_vector
		try:
			initial_index = Line.first_nonzero_index(n)
			terms = [write_coefficient(n[i], is_initial_term=(i==initial_index)) + 'x_{}'.format(i+1) 
					 for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
			output = ' '.join(terms)
		except Exception as e:
			if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
				output = '0'
			else:
				raise e
		constant = round(self.constant_term, num_decimal_places)
		if constant % 1 == 0:
			constant = int(constant)
		return output + ' = {}'.format(constant)
	def __eq__(self, line1):
		if Vector(self.normal_vector).is_zero():
			if not Vector(line1.normal_vector).is_zero():
				return False
			else:
				diff = self.constant_term - line1.constant_term
				return Vector(MyDecimal(diff)).is_near_zero()
		elif Vector(line1.normal_vector).is_zero():
			return False
		if not self.is_parallel(line1):
			return False
		x0 = self.basepoint
		y0 = line1.basepoint
		basepoint_difference = Vector(x0.minus(y0))
		n = Vector(self.normal_vector)
		return basepoint_difference.is_orthogonal_to(n)
	@staticmethod
	def first_nonzero_index(iterable):
		for k, item in enumerate(iterable):
			if not MyDecimal(item).is_near_zero():
				return k
		raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)
	def is_parallel(self, line1):
		u1 = self.normal_vector
		u2 = line1.normal_vector
		u1_vect = Vector(u1)
		u2_vect = Vector(u2)
		return u1_vect.is_parallel_to(u2_vect)
	def intersection(self, line1):
		A = self.normal_vector[0]
		B = self.normal_vector[1]
		C = line1.normal_vector[0]
		D = line1.normal_vector[1]
		k1 = self.constant_term
		k2 = line1.constant_term
		x = (D*k1 - B*k2) / (A*D - B*C)
		y = (-C*k1 + A*k2) / (A*D - B*C)
		try:
			return Vector([x, y]).times_scalar(Decimal('1')/(A*D-B*C))
		except ZeroDivisionError:
			if self == line1:
				return self
			else:
				return None


class MyDecimal(Decimal):
	def is_near_zero(self, eps=1e-10):
		return abs(self) < eps