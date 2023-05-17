from odesolver import ODESolver
from runge_kutta import RungeKuttaMethod


# Define a subclass for Adams method
class AdamsMethod(ODESolver):
    def __init__(self, f, a, b, y0, h, eps, exact_solution):
        super().__init__(f, a, b, h, y0, eps)
        self.exact_solution = exact_solution

    def solve(self):
        # Method to solve the ODE using Adams method and adjust the step size if needed
        n = int((self.b - self.a) / self.h)
        x = self.a
        dots = RungeKuttaMethod(self.f, self.a, self.a + 3 * self.h, self.eps).solve()
        for _ in range(3, n):
            y1 = dots[-1][1] + self.h / 24 * (55 * self.f(dots[-1][0], dots[-1][1]) -
                                              59 * self.f(dots[-2][0], dots[-2][1]) +
                                              37 * self.f(dots[-3][0], dots[-3][1]) -
                                              9 * self.f(dots[-4][0], dots[-4][1]))
            y2 = dots[-1][1] + self.h / 24 * (9 * self.f(x + 4 * self.h, y1) +
                                              19 * self.f(dots[-1][0], dots[-1][1]) -
                                              5 * self.f(dots[-2][0], dots[-2][1]) +
                                              self.f(dots[-3][0], dots[-3][1]))
            x += self.h
            y_exact = [self.exact_solution(x) for x, _ in dots]
            err = max(abs(y_exact[i] - y) for i, (_, y) in enumerate(dots))
            if err > self.eps:
                self.h /= 2
                AdamsMethod(self.f, self.a, self.b, self.y0, self.h, self.epss, self.exact_solution).solve()
            dots.append((x, y2))
        return dots
