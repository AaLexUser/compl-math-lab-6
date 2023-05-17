from data_getter import DataGetter
from euler import EulerMethod
from runge_kutta import RungeKuttaMethod
from adams import AdamsMethod
from get_plot import Plot

if __name__ == '__main__':
    data_getter = DataGetter()
    method, y0, a, b, h, eps, f, exact_f = data_getter.get_data()
    euler = EulerMethod(f=f, a=a, b=b, h=h, y0=y0, eps=eps)
    rk4 = RungeKuttaMethod(f=f, a=a, b=b, h=h, y0=y0, eps=eps)
    adams = AdamsMethod(f=f, a=a, b=b, h=h, y0=y0, eps=eps, exact_solution=exact_f)
    Plot.plot_all(a, b, exact_f, euler.solve(), rk4.solve(), adams.solve())
