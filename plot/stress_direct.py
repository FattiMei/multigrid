import sys
import csvparser
import matplotlib.pyplot as plt


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__} <csv_file_path>")
        sys.exit(1)

    columns, df = csvparser.parse(sys.argv[1])

    plt.figure(1)
    plt.title("Stress test for direct methods")
    plt.xlabel(columns[0])
    plt.ylabel("$||r||$")
    plt.loglog(df[columns[0]], df[columns[1]])
    plt.show()

    plt.figure(2)
    plt.title("Wall time for direct methods")
    plt.xlabel(columns[0])
    plt.ylabel(columns[2])
    plt.loglog(df[columns[0]], df[columns[2]])
    plt.show()
