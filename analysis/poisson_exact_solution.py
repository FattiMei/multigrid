import sys
import sympy as sym
from sympy.utilities.codegen import codegen


# know issue: if the symbolic expression does not contain x or y,
# the generated c function won't have it as an argument
# more research could be done on this issue, for now pay attention to the generated code

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: python {__file__} <src output path> <include output path>')
        sys.exit(1)


    x, y = sym.symbols('x y')

    solution_1d     = sym.exp(-x*x*sym.sin(x / 10))
    forcing_term_1d = -sym.diff(solution_1d,x,2)

    # TODO: fix the issue described above
    solution_2d     = sym.exp(x)*sym.exp(-2*y)
    forcing_term_2d = -sym.diff(solution_2d,x,2) - sym.diff(solution_2d,y,2)

    [(c_name, c_code), (h_name, c_header)] = codegen(
        [
            ('symbolic', 0),
            ('solution_1d', solution_1d),
            ('forcing_term_1d', forcing_term_1d),
            ('solution_2d', solution_2d),
            ('forcing_term_2d', forcing_term_2d)
        ],
        'C99',
        header=True,
        empty=True
    )

    with open(sys.argv[1], 'w') as out:
        out.write(c_code)

    with open(sys.argv[2], 'w') as out:
        out.write(c_header)
