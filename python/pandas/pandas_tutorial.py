# https://realpython.com/pandas-python-explore-dataset/

import numpy as np
import pandas as pd


nba = pd.read_csv('nba_all_elo.csv')
type(nba)

# Also available: read_json(), read_html(), read_sql_table()

len(nba)
nba.shape
nba.head()
nba['gameorder'].head()

pd.set_option('display.precision', 2)
nba.head()
nba.tail(3)

# What are the column data types?
nba.info()
# The "object" data type is a catch all for things that aren't recognized
# This usually means it's a column of strings
# Putting weird things in a pandas df can make it slow

nba.describe()  # This is like summary() from R
nba.describe(include=np.object)  # Includes non-numeric columns
nba['team_id'].value_counts()  # This is like table()
nba['fran_id'].value_counts()
# Find out why Lakers isn't the same count as LAL
nba.loc[nba['fran_id'] == 'Lakers', 'team_id'].value_counts()
# There's another team called the Minneapolis Lakers!!!
nba.loc[nba['team_id'] == 'MNL', 'date_game'].min()
nba.loc[nba['team_id'] == 'MNL', 'date_game'].max()
nba.loc[nba['team_id'] == 'MNL', 'date_game'].agg(('min', 'max'))

# How many points have Boston Celtics scored during all matches?
nba.loc[nba['team_id'] == 'BOS', 'pts'].sum()

# Series objects
revenues = pd.Series([5555, 7000, 1980])
# Series wraps sequence of values + sequences of identifiers (index)
revenues.values
revenues.index
# These are each Numpy 1-D arrays
# The basic difference compared to Numpy is the ability to have arbitrary index values
# i.e., you can use strings as labels in an index
city_revenues = pd.Series(
    [4200, 8000, 6500],
    index=['Amsterdam', 'Toronto', 'Tokyo']
)
city_revenues.values
city_revenues.index
city_revenues
# You can create a series from a dictionary
city_employee_count = pd.Series({'Amsterdam': 5, 'Tokyo': 8})
city_employee_count
city_revenues.keys()
'Tokyo' in city_revenues

# DataFrame objects
city_data = pd.DataFrame({
    'revenue': city_revenues,
    'employee_count': city_employee_count
})
city_data
city_data.shape
# Note that it filled in a missing value with NaN
city_data.index
city_data.values  # This is a Numpy array
city_data.axes

city_revenues[1]
city_revenues['Toronto']
city_revenues[-1]

# .loc gets data via label index
# .iloc gets data via positional index
city_revenues.iloc[1]
city_revenues.loc['Toronto']
# For production code, .loc and .iloc are faster than using regular Python indexing
# .iloc excludes the closing element in a slice (just like Python)
# .loc includes the closing element (weird!)

type(city_data['revenue'])
# You can use the column names as attributes if they're strings
# But this breaks if the column name collides with a method/attribute of the df
city_data.revenue
toys = pd.DataFrame([
    {"name": "ball", "shape": "sphere"},
    {"name": "Rubik's cube", "shape": "cube"}
])
toys["shape"]
toys.shape
# Don't use this for production code!!!
# So I guess I'm just not going to use it, this seems stupid

# For DataFrames, iloc and loc get data row-wise
city_data.iloc[1]
city_data.loc[1]
city_data.loc['Toronto']
city_data.loc['Tokyo': 'Toronto']
# The second argument can select a subset of columns
city_data.loc['Tokyo': 'Toronto', 'revenue']

# Filtering
current_decade = nba[nba['year_id'] > 2010]  # This uses a boolean array to index
current_decade.shape
current_decade

games_with_notes = nba[nba['notes'].notnull()]
games_with_notes.shape
nba[nba['notes'].notna()].shape
# notnull and notna are the same

ers = nba[nba['fran_id'].str.endswith('ers')]
ers
ers['fran_id']

# Baltimore games where both teams scored > 100 points
nba[
    (nba['_iscopy'] == 0) &
    (nba['pts'] > 100) &
    (nba['opp_pts'] > 100) &
    (nba['team_id'] == 'BLB')
]

nba[
    (nba['_iscopy'] == 0) &
    (nba['pts'] > 100) &
    (nba['opp_pts'] > 100) &
    (nba['team_id'] == 'BLB')
]

# Find my Pistons!
nba['fran_id'].value_counts()
pistons = nba.loc[nba['fran_id'] == 'Pistons']
pistons.describe()
pistons['pts'].describe()
pistons['team_id'].value_counts()

# Grouping and aggregating
city_revenues
city_revenues.sum()
city_data.sum()
nba.groupby('fran_id', sort=False)['pts'].sum()
# Losses vs. wins since 2010
nba[
    (nba["fran_id"] == "Pistons") &
    (nba["year_id"] > 2010)
].groupby(["year_id", "game_result"])["game_id"].count()

# Create a copy
df = nba.copy()
df.shape
df['difference'] = df['pts'] - df['opp_pts']
df.shape
df['difference']
renamed_df = df.rename(
    columns={'game_result': 'result',
             'game_location': 'location'}
)
renamed_df.info()
elo_columns = [x for x in df.columns if 'elo' in x]
df.drop(elo_columns, inplace=True, axis=1)
df.shape
df.info()
df['date_game'] = pd.to_datetime(df['date_game'])
df.info()
df['game_location'] = pd.Categorical(df['game_location'])
df.info()

# You can remove rows containing NaN values with .dropna()
nba[nba['pts'] == 0][['notes', 'fran_id', 'opp_fran', 'date_game']]