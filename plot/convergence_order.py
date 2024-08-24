import sys
import csvparser
import matplotlib.pyplot as plt


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__} <csv_file_path>")
        sys.exit(1)

    columns, df = csvparser.parse(sys.argv[1])

    plt.title("Convergence order")
    plt.xlabel(columns[0])
    plt.ylabel("$||err||^{\infty}$")

    for i in range(1,len(columns)):
        plt.loglog(df[columns[0]], df[columns[i]], label=columns[i])

    plt.legend()
    plt.show()
