import json
import re


def get_smoother_abbreviation(smoother):
    abbr = {
        'Jacobi'     : 'JC',
        'GaussSeidel': 'GS',
        'SOR'        : 'SOR',
        'RedBlack'   : 'RB'
    }

    try:
        return abbr[smoother]
    except:
        return smoother


def process_name(complete_name):
    m = re.findall(
        r'.*<(.*?),MgCycle::(.*?),UpdateStrategy::(.*?)(,.*?)?>/(.*)',
        complete_name
    )

    depth, cycle, smoother, threads, n = m[0]
    
    if threads != '':
        threads = int(threads[1:])
        name = f'{cycle}({depth})-{get_smoother_abbreviation(smoother)} ({threads} threads)'
    else:
        name = f'{cycle}({depth})-{get_smoother_abbreviation(smoother)}'

    return (name, int(n))



def parse(json_path):
    result = {}

    with open(json_path, 'r') as file:
        data = json.load(file)

    for bench in data['benchmarks']:
        name, n = process_name(bench['name'])
        time    = bench['real_time']

        try:
            result[name]['n'].append(n)
            result[name]['time'].append(time)

        except KeyError:
            result[name] = {}
            result[name]['n'] = [n]
            result[name]['time'] = [time]

    return result


if __name__ == '__main__':
    print(process_name("BM_step<2,MgCycle::V,UpdateStrategy::GaussSeidel>/4"))
    print(process_name("BM_step<2,MgCycle::V,UpdateStrategy::GaussSeidel,2>/4"))
