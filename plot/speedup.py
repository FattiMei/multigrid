import sys
import jsonparser
import matplotlib.pyplot as plt


def compute_speedup(reference, alternative):
    return [t / s for (s,t) in zip(alternative, reference)]


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {__file__} <csv_file_path>")
        sys.exit(1)

    data = jsonparser.parse(sys.argv[1])

    plt.title('Walltime')
    plt.xlabel('n')
    plt.ylabel('speedup')

    already_printed_reference = False

    for label, perf in data.items():
        if not already_printed_reference:
            reference = perf['time']

            plt.plot(perf['n'], compute_speedup(reference, perf['time']), label=f'{label} (reference)')
            already_printed_reference = True
        else:
            plt.plot(perf['n'], compute_speedup(reference, perf['time']), label=label)


    plt.legend()
    plt.show()
