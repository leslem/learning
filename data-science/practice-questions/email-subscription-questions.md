# Practice questions from email subscriptions

Working through questions from [Interview Query](https://www.interviewquery.com) and [Data Interview Qs](https://www.interviewqs.com) email subscriptions.

Navigate to my Git repo for learning resources and activate my Python `virtualenv`

```
cd ~/devel/learning/data-science/practice-questions
source ~/virtualenv/py-tutorials/bin
```

## Reddit posts and comments
From Data Interview Qs

Suppose you're working for Reddit as an analyst. Reddit is trying to optimize its server allocation per subreddit, and you've been tasked with figuring out how much comment activity happens once a post is published.

Use your intuition to select a timeframe to query the data, as well as how you would want to present this information to the partnering team. The solution will be a SQL query with assumptions that you would need to state if this were asked in an interview. You have the following tables:

Table: posts
| Column Name | Data Type | Description |
| ----------- | --------- | ----------- |
| id | integer | id of the post |
| publisher_id | integer | id the user posting |
| score | integer | score of the post |
| time | integer | post publish time in unix time |
| title | string | title of the post |
| deleted | boolean | is the post deleted? |
| dead | boolean | is the post active? |
| subreddit_id | integer | id of the subreddit |

Table: comments
Column Name	Data Type	Description
id	integer	id of the comment
author_id	integer	id of the commenter
post_id	integer	id of the post the comment is nested under
parent_comment	integer	id of parent comment that comment is nested under
deleted	integer	is comment deleted?

Given the above, write a SQL query to highlight comment activity by subreddit. This problem is intended to test how you can think through vague/open-ended questions. 

### My answer

Main question:
- How much comment activity per ~~post~~ subreddit?

Clarifying questions to ask:
- Comments per post per subreddit? (would normalize by number of posts)
- Comments per subreddit? -- probably more appropriate for this question
- Should I include delete comments or not? Should I include deleted posts or not?
- Should number of unique comment authors be accounted for?
    - Commenting authors per subreddit
    - Commenting authors per post per subreddit
- Why account for timeframe?
    - I think this means only look at recent posts to see the recent comment activity

Approach:
- Join `comments` and `posts` together on `post_id` in order to get `subreddit_id` for each comment
- Total comments per subreddit - aggregation by count
- Include deleted comments, but not deleted posts
- Count all comments by all authors, regardless of uniqueness
- Count comments on posts within the past week

Query:
```
SELECT COUNT(comments.id) AS comment_count, posts.subreddit_id AS subreddit_id
    FROM comments
    LEFT JOIN posts
        ON comments.post_id = posts.id
    WHERE posts.deleted IS NOT TRUE
    GROUP BY subreddit_id
    ORDER BY comment_count DESC
;

# Only include comments within first day after post publication
SELECT COUNT(comments.id) AS comment_count, posts.subreddit_id AS subreddit_id
    FROM comments
    LEFT JOIN posts
        ON comments.post_id = posts.id
    WHERE posts.deleted IS NOT TRUE AND posts.time > NOW() - INTERVAL 1 WEEK
    GROUP BY subreddit_id
    ORDER BY comment_count DESC
;
```

## TV and hypertension case study - BMI
From Data Interview Qs

We'll cover a series of questions over the coming weeks around the following study exploring the relationship between TV consumption and hypertension. Everything you'll need to know about the study is included below, and will be included in future questions.

Background:

Television viewing is strongly associated with an increased risk of childhood and adolescent obesity. However, the association between TV viewing and hypertension in children is unknown. This study aimed to identify whether TV watching is associated with hypertension in obese children.

Methods:

Children seen for obesity, aged 4 to 17 years, were evaluated at three pediatric centers from 2003 to 2005. In 2006–2007, a logistic regression model estimated the odds of hypertension for hours of daily TV time controlling for race, site, and body mass index (BMI) z-score.

Results:

A total of 546 subjects, with a metn age of 12 years, were evaluated. The children had a mean BMI of 35.5±9.3 kg/m2 (98.7th±0.8 percentile, z-score 2.54±0.4). TV time was positively correlated with the severity of obesity. After controlling for race, site, and BMI z-score, both the severity of obesity and daily TV time were significant independent predictors of the presence of hypertension. Children watching 2 to 4 hours of TV had 2.5 times the odds of hypertension compared with children watching 0 to 2 hours. The odds of hypertension for children watching 4 or more hours of TV were 3.3 times greater than for children watching 0 to 2 hours of TV.

Data table:

Table 1 Clinical and demographic characteristics of the study population
Characteristic	All subjects (N=546)	Subjects without hypertension (n=311)	Subjects with hypertension (n=235)
Age, mean (SD), years	11.9 (3.4)	11.7 (3.5)	12.3 (3.3)
Gender, N (%)
 Male	275 (50.4)	160 (51.4)	115 (48.9)
 Female	271 (49.6)	151 (48.6)	120 (51.1)
Race/ethnicity, N (%)
 African American	127 (23.3)	64 (20.6)	63 (26.8)
 Asian/Pacific Islander	20 (3.7)	12 (3.9)	8 (3.4)
 Hispanic	37 (6.8)	23 (7.4)	14 (6.0)
 Multiracial	21 (3.8)	15 (4.8)	6 (2.6)
 White	287 (52.6)	159 (51.1)	128 (54.5)
 Other	54 (9.9)	38 (12.2)	16 (6.8)
Blood pressure (mm Hg)
 Systolic (mean, SD)a	121.0 (16.1)	110.6 (10.4)	134.7 (11.2)
 Diastolic (mean, SD)a	65.3 (9.5)	61.9 (8.4)	69.8 (8.9)
BMI (kg/m2)
 Mean (SD)a	35.5 (9.3)	33.9 (8.5)	37.6 (9.8)
 Percentile, mean (SD)a	98.7 (0.8)	98.6 (0.9)	98.9 (0.5)
 Z-score, mean (SD)a	2.54 (0.4)	2.49 (0.5)	2.59 (0.4)
Hours of TV/day, Mean (SD)a	3.1 (1.8)	2.8 (1.7)	3.6 (1.8)
TV time categories, N (%)
 0 to <2 hours/day	121 (22.2)	94 (30.2)	27 (11.5)
 2 to <4 hours/day	202 (37.0)	115 (37.0)	87 (37.0)
 ≥4 hours/day	223 (40.8)	102 (32.8)	121 (51.5)

Question: Refer to the data table shown above, and develop a hypothesis test to compare the mean BMI of subjects without hypertension to the mean BMI of subjects with hypertension (you're testing to see if there's a significant difference in the means). Walk through what test you chose and why, run your test (feel free to use tech here, doing it by hand is not needed), and share your conclusion.

### My answer

Clarifying questions:
- So we want to know if TV watching is associated with hypertension in obsese children? -- no, they already did this
- What they did was a logsitic regression
- I think these are the two models they evaluated:
odds(hypertension) ~ hours of TV time + race + site + BMI z-score
odds(hypertension) ~ obesity severity + race + site + BMI z-score
- we want to do something simpler - just comparing the mean BMI in subjects with vs. without hypertension

My approach:
- t-test, 2-tailed

```{r}
bmi_test <-
t.test(x = bmi_in_hypertension,
       y = bmi_in_nonhypertension,
       alternative=('two.sided'))
```

Output would look like:
bmi_test$statistic
bmi_test$p.value
bmi_test$conf.int

If p.value is ≤ 0.05, then reject the null hypothesis of equal means between the two groups. If p.value is > 0.05, then there is no statistically significant difference in BMI between the two groups.

## Relationship between fitness and smoking
From Data Interview Qs

Given the data table below, determine if there is a relationship between fitness level and smoking habits:

Low fitness level	Medium-low fitness level	Medium-high fitness level	High fitness level
Never smoked	113	113	110	159
Former smokers	119	135	172	190
1 to 9 cigarettes daily	77	91	86	65
>=10 cigarettes daily	181	152	124	73

You don't have to fully solve for the number here (that would be pretty time-intensive for an interview setting), but lay out the steps you would take to solve such a problem. 

### My answer
I'm being asked to determine if there's a relationship between frequencies of two categorical variables.
I'm given a contingency table.
I need to use this table as input to a Chi squared contingency table test.
"Are my observed counts in the cells of the table what I would expect by chance if there's no relationship between the two categorical variables."
Degrees of freedom 
