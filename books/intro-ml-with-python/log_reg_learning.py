# Example from near Figure 2-15 of Introduction to Machine Learning with Python

import maptlotlib.pyplot as plt

from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix, accuracy_score

skcancer = load_breast_cancer()
# Split data into testing and training, stratifying on the target variable to ensure balance
X_train, X_test, y_train, y_test = train_test_split(
    skcancer.data, skcancer.target, stratify=skcancer.target, random_state=42)

classifier = LogisticRegression().fit(X_train, y_train)
# This gives me a "failed to converge" warning

# score() gives mean accuracy on a given dataset
print('Training set score: {:.3f}'.format(classifier.score(X_train, y_train)))  # 0.946
print('Test set score: {:.3f}'.format(classifier.score(X_test, y_test)))  # 0.958

plt.plot(classifier.coef_.T, 'o')
plt.show()

# Look at coefficients
dir(skcancer)
skcancer.feature_names
classifier.coef_


# K-fold cross-validation
from sklearn.model_selection import cross_val_score

scores = cross_val_score(LogisticRegression(), skcancer.data, skcancer.target, cv=10)
scores

# Grid search to tune parameters
from sklearn.model_selection import GridSearchCV

param_grid = {
    'C': 
}

grid_search = GridSearchCV(LogisticRegression(), param_grid, cv=5)


# Not sure how to get at the parameter values that were tried...
classifier = LogisticRegressionCV(cv=5, Cs=10).fit(skcancer.data, skcancer.target)
dir(classifier)
classifier.scores_
classifier.get_params()
classifier._get_param_names()
