# Data science practice problems from Jay Feng
[Notes and technical questions from interviewing as a Data Scientist in 2018](https://towardsdatascience.com/notes-and-technical-questions-from-interviewing-as-a-data-scientist-in-2018-20e7e3ee4ab3)

## Coding
* Fizzbuzz
    * http://wiki.c2.com/?FizzBuzzTest
    * Write a program that prints the numbers from 1 to 100. But for multiples of three print “Fizz” instead of the number and for the multiples of five print “Buzz”. For numbers which are multiples of both three and five print “FizzBuzz”.

```{r}
# Multiples of 3: replace with "fizz"
# Multiples of 5: replace with "buzz"
# Multiples of both: replace with "fizzbuzz"
for (x in seq(100)) {
    out <- ""
    if (x %% 3 == 0) {
        out <- paste0(out, "Fizz")
    }
    if (x %% 5 == 0) {
        out <- paste0(out, "Buzz")
    }
    if (out == "") {
        out <- x
    }
    # print(paste(x, ":", out))
    print(out)
}
```

I think this is a more idiomatic way to write this in R
```{r}
df <- data.frame(x=seq(100))
df$div3 <- df$x %% 3 == 0
df$div5 <- df$x %% 5 == 0
df$out <- ifelse(df$div3 & df$div5, 'FizzBuzz', ifelse(df$div3, 'Fizz', ifelse(df$div5, 'Buzz', '')))
```

```{python}
for x in range(1, 100):
    out = ''
    if x % 3 == 0:
        out += 'Fizz'
    if x % 5 == 0:
        out += 'Buzz'
    if out == '':
        out = str(x)
    print(out)
```

* Given a list of timestamps in sequential order, return a list of lists grouped by weekly aggregation.

```{r}
library(lubridate)
library(tidyverse)
library(data.table)

timestamps <- sample(seq(as.POSIXct('2018-01-01'), as.POSIXct('2018-12-31'), 1), 100)
timestamps <- sort(timestamps)
week_nums <- week(timestamps)
df <- data.frame(timestamp=timestamps, week_num=week_nums)
# You can use this dplyr function to make the list of lists by week_num.
# But I really think it's easiest to keep it as a dataframe with week_num labels.
df %>% nest(timestamp)

# For example, it's easier just to use data.table to aggregate things by week!
dt <- data.table(df)
dt
dt[ , .(.N), by=week_num]
```

```{python}
from datetime import datetime, timedelta
import random
from collections import defaultdict

start = datetime(2018, 1, 1, 00, 00, 00)
end = datetime(2018, 12, 31, 00, 00, 00)
year_delta = end - start
# You can't multiply a datetime by a fraction (from random.random)
# but you can multiply a timedelta by a fraction, so just add a random fraction
# of the timedelta for one year to the starting datetime.
# Source: https://gist.github.com/rg3915/db907d7455a4949dbe69
timestamps = [start + (year_delta * random.random()) for x in range(100)]
week_timestamps = [(x.isocalendar()[1], x) for x in timestamps]

# Use the defaultdict from collections to make a dict with list as the default factory
# This appends all of the values with the same week number together into a list
timestamps_by_week = defaultdict(list)
for k, v in week_timestamps:
    timestamps_by_week[k].append(v)

# If you really want the timestamps collected together into lists by week,
# then I think this is a better data structure than a plain list of lists.
# But just as in the R example, it may be better to just leave them as
# columns in a data frame (e.g. pandas) and aggregate on the week column.
```

* Given a list of characters, a list of prior of probabilities for each character, and a matrix of probabilities for each character combination, return the optimal sequence for the highest probability.

```{r}
# Set up a simple example dataset.
chars <- c('a', 'b', 'c')
probs <- runif(length(chars))
probs <- probs / sum(probs)
combo_probs <- matrix(runif(length(chars) ** 2), nrow=3)
combo_probs <- combo_probs / sum(combo_probs)
rownames(combo_probs) <- chars
colnames(combo_probs) <- chars

```


* Given a log file with rows featuring a date, a number, and then a string of names, parse the log file and return the count of unique names aggregated by month.
## Product
* Given there are no metrics being tracked for Google Docs, a product manager comes to you and asks what are the top five metrics you would implement?
* In addition, let’s say theres a dip in the engagement metric of Google Docs. What would you investigate?
* Let’s say we want to implement a notification system for reminding nurses to discharge patients at a hospital. How would you implement it?
* Let’s say at LinkedIn we want to implement a green dot for an “active user” on the new messaging platform. How would you analyze the effectiveness of it for roll out?
## SQL
* Given a payment transactions table and a customers table, return the customer’s name and the first transaction that the customer made.
* Given a payments transactions table, return a frequency distribution of the number of payments each customer made. (I.E. 1 transaction — 100 customers, 2 transactions — 50 customers, etc…)
* Given the same payments table, return the cumulative distribution. (At least one transaction, at least two transactions, etc…)
* Given a table of — friend1|friend2. Return the number of mutual friends between two friends.
## AB Testing
* Given AB test funnel statistics such as the sample size, sign up rate, feature 1 usage rate, feature 2 usage rate, analyze which variant won and why.
* How would you design an experiment to change a button on a sign up page?
* How do you know if you have enough sample size?
* How do you run significance tests on more than one variant?
* How do you reduce variance and bias in an AB test?
* Explain a P-value and confidence interval to a product manager or non-technical person.
## Machine Learning
* What features would you use to predict the time spent for a restaurant preparing food from the moment an order comes in?
* Can you come up with a scenario in which you would rather under-predict versus over-predict?
* Analyzing the results of a model, how would you explain the tradeoff between bias and variance?
* Explain how a Random Forest model actual works under the hood.
* How do you know if you have enough data for your model?
* How do you evaluate a model? (F1 score, ROC curve, cross validation, etc…)
## Probability
* Given uniform distributions X and Y and the mean 0 and standard deviation 1 for both, what’s the probability of 2X > Y?
* There are four people in an elevator and four floors in a building. What’s the probability that each person gets off on a different floor?
* What’s the probability that two people get off on the same floor?
* Given a deck of cards labeled from 1 to 100, what’s the probability of getting Pick 1 < Pick2 < Pick3?