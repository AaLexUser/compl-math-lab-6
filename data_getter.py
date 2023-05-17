import math
class DataGetter:
    def __init__(self):
        # Initialize the attributes of the class
        self.method = None  # ODE method
        self.y0 = None  # initial value of y
        self.a = None  # lower bound of x
        self.b = None  # upper bound of x
        self.h = None  # step size
        self.eps = None  # error tolerance
        self.f = None  # function defining the ODE
        self.exact_f = None  # exact solution of the ODE

    def get_method(self):
        # Method to get the ODE method from the user
        print('Choose ODE method:')
        print('1. Euler method')
        print('2. Runge-Kutta method')
        print('3. Adams method')
        self.method = int(input())
        while self.method < 1 or self.method > 3:
            print('Wrong input. Try again:')
            self.method = int(input())

    def get_initial_value(self):
        # Method to get the initial value of y from the user
        print('Enter y0:')
        self.y0 = float(input())

    def get_interval(self):
        # Method to get the interval of x from the user
        print('Enter a, b:')
        while True:
            self.a, self.b = map(float, input().split())
            if self.a >= self.b:
                print('Wrong input. Try again:')
            else:
                break

    def get_step_size(self):
        # Method to get the step size from the user
        print('Enter h:')
        while True:
            self.h = float(input())
            if self.h <= 0 or self.h > self.b - self.a:
                print('Wrong input. Try again:')
            else:
                break

    def get_error_tolerance(self):
        # Method to get the error tolerance from the user
        print('Enter eps:')
        while True:
            self.eps = float(input())
            if self.eps <= 0:
                print('Wrong input. Try again:')
            else:
                break

    def get_function(self):
        # Method to get the function defining the ODE from the user
        print('Choose function:')
        print('1. y\' = y - x^2 + 1  на [1, 2]  y(1) = 2')
        print('2. y\' = 1 + (y - x)^2 на [0, 1]  y(0) = 1')
        print('3. y\' = y / x + x / y на [1, 2]  y(1) = 1')
        print('4. y\' = y + (1 + x) * y^2 на [1, 1.5]  y(1) = -1')
        functions = [lambda x, y: y - x ** 2 + 1,
                     lambda x, y: 1 + (y - x) ** 2,
                     lambda x, y: y / x + x / y,
                     lambda x, y: y + (1 + x)* y ** 2]
        exact_f = [
            lambda x: (x + 1) ** 2 - 2 * math.exp(x - 1),
            lambda x: (x**2 - x - 1) / (x - 1),
            lambda x: x * math.sqrt(2 * math.log(x) + 1),
            lambda x: -1 / (x),


        ]
        self.exact_f = exact_f[int(input()) - 1]
        self.f = functions[int(input()) - 1]

    def get_data(self):
        # Method to get all the data from the user by calling other methods
        self.get_method()
        self.get_initial_value()
        self.get_interval()
        self.get_step_size()
        self.get_error_tolerance()
        self.get_function()
        return self.method, self.y0, self.a, self.b, self.h, self.eps, self.f, self.exact_f
