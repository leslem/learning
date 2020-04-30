# Questions from Data Science Central post

[The original post](https://www.datasciencecentral.com/profiles/blog/show?id=6448529:BlogPost:371610)

## Questions from Vincent Granville
- What is the life cycle of a data science project?
- How do you measure yield (over base line) resulting from a new or refined algorithm or architecture?
- What is cross-validation? How to do it right?
- Is it better to design robust or accurate algorithms?
- Have you written production code? Prototyped an algorithm? Created a proof of concept?
- What is the biggest data set you have worked with, in terms of training set size, and in terms of having your algorithm implemented in production mode to process billions of transactions per day / month / year?
- Name a few famous API's (for instance Google search). How would you create one?
- How to efficiently scrape web data, or collect tons of tweets?
- How to optimize algorithms (parallel processing and/or faster algorithm: provide examples for both)
- Examples of NoSQL architecture?
- How do you clean data?
- How do you define / select metrics? Have you designed and used compound metrics?
- Examples of bad and good visualizations?
- Have you been involved - as an adviser or architect - in the design of dashboard or alarm systems?
- How frequently an algorithm must be updated? What about lookup tables in real-time systems?
- Provide examples of machine-to-machine communication.
- Provide examples where you automated a repetitive analytical task.
- How do you assess the statistical significance of an insight?
- How to turn unstructured data into structured data?
- How to very efficiently cluster 100 billion web pages, for instance with a tagging or indexing algorithm? 
- If you were interviewing a data scientist, what questions would you ask her?

## Questions from Jay Verkuilen on Quora
- Explain what regularization is and why it is useful. What are the benefits and drawbacks of specific methods, such as ridge regression and LASSO?
- Explain what a local optimum is and why it is important in a specific context, such as k-means clustering. What are specific ways for determining if you have a local optimum problem? What can be done to avoid local optima?
- Assume you need to generate a predictive model of a quantitative outcome variable using multiple regression. Explain how you intend to validate this model.
- Explain what precision and recall are. How do they relate to the ROC curve?
- Explain what a long tailed distribution is and provide three examples of relevant phenomena that have long tails. Why are they important in classification and prediction problems?
- What is latent semantic indexing? What is it used for? What are the specific limitations of the method?
- What is the Central Limit Theorem? Explain it. Why is it important? When does it fail to hold?
- What is statistical power?
- Explain what resampling methods are and why they are useful. Also explain their limitations.
- Explain the differences between artificial neural networks with softmax activation, logistic regression, and the maximum entropy classifier.
- Explain selection bias (with regards to a dataset, not variable selection). Why is it important? How can data management procedures such as missing data handling make it worse?
- Provide a simple example of how an experimental design can help answer a question about behavior. For instance, explain how an experimental design can be used to optimize a web page. How does experimental data contrast with observational data.
- Explain the difference between "long" and "wide" format data. Why would you use one or the other?
- Is mean imputation of missing data acceptable practice? Why or why not?
- Explain Edward Tufte's concept of "chart junk." 
- What is an outlier? Explain how you might screen for outliers and what you would do if you found them in your dataset. Also, explain what an inlier is and how you might screen for them and what you would do if you found them in your dataset.
- What is principal components analysis (PCA)? Explain the sorts of problems you would use PCA for. Also explain its limitations as a method.
- You have data on the duration of calls to a call center. Generate a plan for how you would code and analyze these data. Explain a plausible scenario for what the distribution of these durations might look like. How could you test (even graphically) whether your expectations are borne out?
- Explain what a false positive and a false negative are. Why is it important to differentiate these from each other? Provide examples of situations where (1) false positives are more important than false negatives, (2) false negatives are more important than false positives, and (3) these two types of errors are about equally important.
- Explain likely differences encountered between administrative datasets and datasets gathered from experimental studies. What are likely problems encountered with administrative data? How do experimental methods help alleviate these problems? What problems do they bring?

## Questions from Kavita Ganesan
- What is a gold standard ? 
	- Believe it or not there are data scientists (even at very senior levels) who claim to know a hell lot about supervised machine learning and know nothing about what a gold standard is!
- What is the difference between supervised learning and unsupervised learning? - Give concrete examples.
- What does NLP stand for?
	- Some data scientists claim to also do NLP.  
- Write code to count the number of words in a document using any programming language. Now, extend this for bi-grams.
	- I have seen a senior level data scientist who actually struggled to implement this. 
- What are feature vectors?
- When would you use SVMs vs Random Forrest and Why?
- What is your definition of Big Data, and what is the largest size of data you have worked with? Did you parallelize your code?
	- If their notion of big data is just volume - you may have a problem. Big Data is more than just volume of data. If the largest size of data they have worked with is 5MB - again you may have a problem.
- How do you work with large data sets?
	- If the answer only comes out as hadoop it clearly shows that their view of solving problems is extremely narrow. Large data problems can be solved with:
	1. efficient algorithms
	2. multi-threaded applications
	3. distributed programming
	4. more...
- Write a mapper function to count word frequencies (even if its just pseudo code)
- Write a reducer function for counting word frequencies (even if its just pseudo code)
