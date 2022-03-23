import concurrent.futures


def func(x, y, z):
    return (x * y) + (x * z) + (y * z)

data = [((x, y, z), func(x, y, z)) for x in range(10) for y in range(10) for z in range(10)]

with concurrent.futures.ProcessPoolExecutor(max_workers=24) as executor:
    results = [(args, executor.submit(func,*args).result()) for args, _ in data]
    
print(data == results)
