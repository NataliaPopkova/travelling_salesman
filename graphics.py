import matplotlib.pyplot as plt

with open('Seq.txt') as f:
    seq = [int(x) for x in f]

with open('Cities.txt') as f:
    cities = []
    for line in f:
        cities.append([float(x) for x in line.split()])

cities_seq = []
for number in seq:
    cities_seq.append(cities[number])

x = [element[0] for element in cities_seq]
y = [element[1] for element in cities_seq]


plt.plot(x, y, 'xb-')

plt.show()
