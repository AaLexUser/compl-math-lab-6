import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

"""
    Approximates the solution to a first-order ODE using Euler's method.
    Args:
        f: function defining the ODE dy/dx = f(x, y)
        a, b: interval for x
        h: step size
        y0: initial value of y
        eps: maximum error tolerance
    Returns:
        List of tuples representing the (x, y) coordinates of the approximate solution.
"""
def euler_method(f, a, b, h, y0, eps):
    n = int((b - a) / h)
    x = a
    p = 1
    dots = [(a, y0)]
    for _ in range(n):
        y0 += h * f(x, y0)
        y_h = y0 + h * f(x, y0)
        y_h2 = y0 + (h / 2) * f(x, y0)
        y_h2 = y_h2 + (h / 2) * f(x + h / 2, y_h2)
        err = abs(y_h - y_h2) / (2 ** p - 1)
        while err > eps:
            h /= 2
            y_h = y0 + h * f(x, y0)
            y_h2 = y0 + (h / 2) * f(x, y0)
            y_h2 = y_h2 + (h / 2) * f(x + h / 2, y_h2)
            err = abs(y_h - y_h2) / (2 ** p - 1)
        y0 = y_h2
        x += h
        dots.append((x, y0))
    return dots


def runge_kutta_method(f, a, b, h, y0, eps):
    n = int((b - a) / h)
    x = a
    dots = [(a, y0)]
    for _ in range(n):
        k1 = h * f(x, y0)
        k2 = h * f(x + h / 2, y0 + k1 / 2)
        k3 = h * f(x + h / 2, y0 + k2 / 2)
        k4 = h * f(x + h, y0 + k3)
        y0 += (k1 + 2 * k2 + 2 * k3 + k4) / 6
        x += h
        dots.append((x, y0))
        while err > eps:
            h /= 2
            k1 = h / 2 * f(x, y0)
            k2 = h / 2 * f(x + h / 4, y0 + k1 / 2)
            k3 = h / 2 * f(x + h / 4, y0 + k2 / 2)
            k4 = h * f(x + h / 2, y0 + k3)
            yh2 = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            err = abs(yh - yh2) / (2 ** 4 - 1)
        y0 = yh2
        x += h
        dots.append((x, y0))
    return dots


def adams_method(f, a, b, y0, h, eps, exact_solution):
    n = int((b - a) / h)
    x = a
    dots = runge_kutta_method(f, a, a + 3 * h, h, y0)
    for _ in range(3, n):
        y1 = dots[-1][1] + h / 24 * (55 * f(dots[-1][0], dots[-1][1]) -
                                     59 * f(dots[-2][0], dots[-2][1]) +
                                     37 * f(dots[-3][0], dots[-3][1]) -
                                     9 * f(dots[-4][0], dots[-4][1]))
        y2 = dots[-1][1] + h / 24 * (9 * f(x + 4 * h, y1) +
                                     19 * f(dots[-1][0], dots[-1][1]) -
                                     5 * f(dots[-2][0], dots[-2][1]) +
                                     f(dots[-3][0], dots[-3][1]))
        x += h
        y_exact = [exact_solution(x) for x, _ in dots]
        err = max(abs(y_exact[i] - y) for i, (_, y) in enumerate(dots))
        if err > eps:
            h /= 2
            adams_method(f, a, b, y0, h, eps, exact_solution)
        dots.append((x, y2))
    return dots


def get_data():
    print('Choose ODE method:')
    print('1. Euler method')
    print('2. Runge-Kutta method')
    print('3. Adams method')
    method = int(input())
    while method < 1 or method > 3:
        print('Wrong input. Try again:')
        method = int(input())
    print('Enter y0:')
    y0 = float(input())
    print('Enter a, b:')
    while True:
        a, b = map(float, input().split())
        if a >= b:
            print('Wrong input. Try again:')
        else:
            break
    print('Enter h:')
    while True:
        h = float(input())
        if h <= 0 or h > b - a:
            print('Wrong input. Try again:')
        else:
            break
    print('Enter eps:')
    while True:
        eps = float(input())
        if eps <= 0:
            print('Wrong input. Try again:')
        else:
            break
    print('Choose function:')
    print('1. y\' = y - x^2 + 1')
    print('2. y\' = 1 + (y - x)^2')
    print('3. y\' = y / x + x / y')
    print('4. y\' = y^2 - x^2 + 1')
    functions = [lambda x, y: y - x ** 2 + 1,
                    lambda x, y: 1 + (y - x) ** 2,
                    lambda x, y: y / x + x / y,
            lambda x, y: y ** 2 - x ** 2 + 1]
    function = functions[int(input()) - 1]
    return method, y0, a, b, h, eps, function


if __name__ == '__main__':
    pass
