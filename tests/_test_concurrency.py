import concurrent.futures

def func(args):
    x, y, z = args
    return (x * y) + (x * z) + (y * z)

data = [(x, y, z) for x in range(10) for y in range(10) for z in range(10)]
expected = [(x * y) + (x * z) + (y * z) for x, y, z in data]

num_lines = 12
chunk = []
results = []
for i, args in enumerate(data):
    chunk.append(args)
    if ((i + 1) % num_lines) == 0:
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results_ = executor.map(func, chunk)
        results.extend(results_)
        chunk = []
if chunk:
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results_ = executor.map(func, chunk)
    results.extend(results_)
    chunk = []

failed = False
for r, e in zip(results, expected):
    if r != e:
        failed = True
        print(f'{r} != {e}')

if failed:
    print('Failed')
else:
    print('Passed')