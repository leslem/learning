# # Chapter 1 exercises from *Introduction to Machine Learning with Python* by Müller and Guido
# Exploring the basics of how scikit-learn works.

# +
import matplotlib.pyplot as plt
import mglearn
import numpy as np
import pandas as pd

from sklearn.datasets import load_iris
from sklearn import model_selection
from sklearn import neighbors
# -

# ## Load the iris dataset.

# +
iris = load_iris()
dir(iris)
print(iris['DESCR'])
print(iris['data'])
print(iris['feature_names'])
print(iris['filename'])
print(iris['target'])
print(iris['target_names'])
# -

# ## Split data into test and training sets (25% test)

# +
X_train, X_test, y_train, y_test = model_selection.train_test_split(iris['data'], iris['target'], random_state=0)
print(X_train.shape)
print(y_train.shape)
print(X_test.shape)
print(y_test.shape)
# -

# ## Plot all variables from training data in a pair plot

# +
iris_df = pd.DataFrame(X_train, columns=iris.feature_names)
pd.plotting.scatter_matrix(iris_df, c=y_train, figsize=(15, 15), marker='o',
                           hist_kwds={'bins': 20}, s=60, alpha=0.8,
                           cmap=mglearn.cm3)  # Use custom colors from the book's package
plt.show()
# -

# The Pandas pair plot shows a matrix of feature scatter plots, in this case colored by iris species.
# The diagonals show histograms of each feature rather than plotting it against itself.

# ## Build a nearest neighbors model

# +
knn = neighbors.KNeighborsClassifier(n_neighbors=1)
knn.fit(X_train, y_train)
dir(knn)
# -

# ## Make predictions for unseen data

# +
X_new = np.array([[5, 2.9, 1, 0.2]])
print(X_new.shape)
prediction = knn.predict(X_new)
print('Prediction:', prediction)
print('Predicted target name:', iris['target_names'][prediction])
# You need to map the prediction value to the target names in the weirdly formatted dataset.
# -

# ## Evaluate the model with the test set

# +
y_test_predictions = knn.predict(X_test)
# Accuracy is the percentage of these predictions that match the true labels.
print('Accuracy: {:.2f}'.format(np.mean(y_test_predictions == y_test)))
# But the classifier object has a built in scoring method:
knn.score(X_test, y_test)
# This is compared to the baseline of 0.333, which is what you'd get if you pick one to predict every time.
# -

# ## Summary
# The basic scikit-learn process is as follows:
# 1. Instantiate a model with your desired parameters
# 2. Fit the model to your input data (training set)
# 3. Make predictions with you model (if applicable) on the test set
# 4. Score the performance of your model by comparing predictions to truth on the test set
