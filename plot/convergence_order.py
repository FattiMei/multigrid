import sys
import numpy as np
import csvparser
import matplotlib.pyplot as plt


SUPPOSED_ORDER = -2


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__} <csv_file_path>")
        sys.exit(1)

    columns, df = csvparser.parse(sys.argv[1])

    plt.title("Convergence order")
    plt.xlabel(columns[0])
    plt.ylabel("$||err||^{\infty}$")

    n = np.array(df[columns[0]], dtype=np.float64)

    for i in range(1,len(columns)):
        plt.loglog(n, df[columns[i]], label=columns[i])

    plt.loglog(n, n ** SUPPOSED_ORDER, label=f'O(N^{SUPPOSED_ORDER})')

    plt.legend()
    plt.show()
