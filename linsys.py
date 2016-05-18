from decimal import Decimal, getcontext
from copy import deepcopy
from vector import Vector
from plane import Plane

getcontext().prec = 30

class LinearSystem(object):
	ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
	NO_SOLUTIONS_MSG = 'No solutions'
	INF_SOLUTIONS_MSG = 'Infinitely many solutions'
	def __init__(self, planes):
		try:
			d = planes[0].dimension
			for p in planes:
				assert p.dimension == d
			self.planes = planes
			self.dimension = d
		except AssertionError:
			raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)
	def __str__(self):
		ret = 'Linear System:\n'
		temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
		ret += '\n'.join(temp)
		return ret
	def __len__(self):
		return len(self.planes)
	def __getitem__(self, i):
		return self.planes[i]
	def __setitem__(self, i, x):
		try:
			assert x.dimension == self.dimension
			self.planes[i] = x
		except AssertionError:
			raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)
	def indices_of_first_nonzero_terms_in_each_row(self):
		num_equations = len(self)
		num_variables = self.dimension
		indices = [-1] * num_equations
		for i, p in enumerate(self.planes):
			try:
				indices[i] = p.first_nonzero_index(p.normal_vector)
			except Exception as e:
				if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
					continue
				else:
					raise e
		return indices
	# 1 swapping two rows
	def swap_rows(self, row1, row2):
		self[row1], self[row2] = self[row2], self[row1] 
	# 2 coefficient * row
	def multiply_coefficient_and_row(self, coefficient, row):
		mult = Vector(self[row]).scalar_multiply
		const = coefficient * self.planes[row].constant_term
		self.planes[row] = Plane(mult, const)
	# 3 (coefficient * row1) + row2
	def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
		row_to_add_copy = self.planes[row_to_add] # make copy of row0
		multiplied = [coefficient*x for x in (row_to_add_copy.normal_vector)] # multiply row
		k = coefficient * row_to_add_copy.constant_term # multiply row0's constant
		k1 = self.planes[row_to_be_added_to].constant_term
		instance_row = self.planes[row_to_be_added_to] # make a copy of row1
		new_row_to_be_added_to = [x+y for x,y in zip(multiplied, instance_row.normal_vector)] # row0 + row1
		self.planes[row_to_be_added_to] = Plane(Vector(new_row_to_be_added_to), k + k1) # change row1
	def compute_triangular_form(self):
		system = deepcopy(self)
		num_equations = len(system)
		num_variables = system.dimension
		j = 0
		for i in range(num_equations):
			while j < num_variables:
				c = MyDecimal(system[i].normal_vector[j])
				if c.is_near_zero():
					swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient_if_able(i, j)
					if not swap_succeeded:
						j += 1
						continue
					system.clear_coefficients_below(i, j)
					j += 1
				break
		return system
	def swap_with_row_below_for_nonzero_coefficient_if_able(self, row, col):
		num_equations = len(self)
		for k in range(row+1, num_equations):
			coefficient = MyDecimal(self[k].normal_vector[col])
			if not coefficient.is_near_zero():
				self.swap_rows(row, k)
				return True
		return False
	def clear_coefficients_below(self, row, col):
		num_equations = len(self)
		beta = MyDecimal(self[row].normal_vector[col])
		for k in range(row+1, num_equations):
			n = self[k].normal_vector
			gamma = n[col]
			alpha = -gamma / beta
			self.add_multiple_times_row_to_row(alpha, row, k)
	def compute_rref(self):
		tf = self.compute_triangular_form()
		num_equations = len(tf)
		pivot_indices = tf.indices_of_first_nonzero_terms_in_each_row()
		for i in range(num_equations)[::-1]:
			j = pivot_indices[i]
			if j < 0:
				continue
			tf.scale_row_to_make_coefficient_equal_one(i, j)
			tf.clear_coefficients_above(i, j)
		return tf
	def scale_row_to_make_coefficient_equal_one(self, row, col):
		n = self[row].normal_vector
		beta = Decimal('1.0') / n[col]
		self.multiply_coefficient_and_row(beta, row)
	def clear_coefficients_above(self, row, col):
		for k in range(row)[::-1]:
			n = self[k].normal_vector
			alpha = -(n[col])
			self.add_multiple_times_row_to_row(alpha, row, k)
	def gaussian_elimination_to_extract_solution(self):
		rref = self.compute_rref()
		



class MyDecimal(Decimal):
	def is_near_zero(self, eps=1e-10):
		return abs(self) < eps


a = Plane(Vector(['5.262','2.739','-9.878']), '-3.441')
b = Plane(Vector(['5.111','6.358','7.638']), '-2.152')
c = Plane(Vector(['2.016','-9.924','-1.367']), '-9.278')
d = Plane(Vector(['2.167','-13.543','-18.883']), '-10.567')
s = LinearSystem([a,b,c,d])