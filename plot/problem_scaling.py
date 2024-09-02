import sys
import jsonparser
import matplotlib.pyplot as plt


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__} <csv_file_path>")
        sys.exit(1)

    data = jsonparser.parse(sys.argv[1])

    plt.title('Walltime')
    plt.xlabel('n')
    plt.ylabel('walltime[ns]')

    already_printed_reference = False

    for label, perf in data.items():
        if not already_printed_reference:
            plt.loglog(perf['n'], [n**2 for n in perf['n']], label='O(n^2)')
            already_printed_reference = True

        plt.loglog(perf['n'], perf['time'], label=label)


    plt.legend()
    plt.show()
