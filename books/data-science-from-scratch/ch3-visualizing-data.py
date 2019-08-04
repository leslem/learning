"""Visualizing data."""

from matplotlib import pyplot as plt
import random
from collections import Counter
from string import ascii_lowercase

# Line chart example from Fig 3-1.
years = range(1950, 2020, 10)
gdp = [random.random() * 18000 for _ in years]
gdp.sort()
plt.plot(years, gdp, color='blue', marker='o', linestyle='solid')
plt.title("Nominal GDP")
plt.ylabel("Billions of $")
plt.show()
plt.close()

# Bar chart example from Fig 3-2.
cats = ['Emma', 'Wallace', 'Darwin', 'Oliver']
cuddles_per_day = [3, 25, 18, 8]
plt.bar(range(len(cats)), cuddles_per_day)
plt.title('Cuddles required by each cat')
plt.ylabel('Cuddles per day')
plt.xticks(range(len(cats)), cats)
plt.show()
plt.close()

# Histogram example from Fig 3-3.
grades = [random.choice(range(25, 100)) for _ in range(100)]
histogram = Counter(min(grade // 10 * 10, 90) for grade in grades)
plt.bar([x + 5 for x in histogram.keys()],
        histogram.values(),
        10,
        edgecolor=(0, 0, 0)  # RGB black
        )
plt.axis([-5, 105, 0, round(max(histogram.values())) + 5])
plt.xticks([10 * i for i in range(11)])
plt.xlabel("Grade decile")
plt.ylabel("Number of students")
plt.title("Distribution of exam grades")
plt.show()
plt.close()

# Scatterplot example from Fig 3-7.
friends = [round(random.random() * 500) for _ in range(50)]
minutes = [round(random.random() * 300) for _ in range(len(friends))]
# long_alphabet = ascii_lowercase * (round(len(friends) / len(ascii_lowercase)) + 1)
# labels = [long_alphabet[_] for _ in range(len(friends))]
labels = ["person_{:2d}".format(_) for _ in range(len(friends))]
plt.scatter(friends, minutes)
for label, friend_count, minute_count in zip(labels, friends, minutes):
    plt.annotate(
        label,
        xy=(friend_count, minute_count),
        xytext=(5, -5),
        textcoords='offset points'
    )
plt.title("Daily minutes vs. Number of Friends")
plt.xlabel('Number of friends')
plt.ylabel('Daily minutes spent on the site')
plt.show()
plt.close()
