# # set up virtualenv that has the latest python, numpy, and scipy
# /projects/resources/software/python/bin/virtualenv ~/virtualenv/scipy_py37
# # When I set this virtualenv up in my workon_home directory, it shows up in the listing
# # from virtualenvwrapper
# lsvirtualenv
# # and it works with virtualenvwrapper and still uses the expected version of python
# workon scipy_py37
# which python
# python --version
# # Install scipy etc.
# pip --version
# pip install numpy scipy matplotlib ipython jupyter pandas sympy nose

# # Start the python interpreter
# python

# These are standard import statements and the aliases appear to be community-accepted
# conventions. I'm not sure how I feel about them in production code though.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

a = np.array(range(5))

# Compare to timing of python data structures
L = list(range(1000))
%timeit [i**2 for i in L]
# 449 µs ± 193 ns per loop (mean ± std. dev. of 7 runs, 1000 loops each)
a = np.arange(1000)
%timeit a**2
# 4.02 µs ± 41.1 ns per loop (mean ± std. dev. of 7 runs, 100000 loops each)

# Two orders of magnitude speed increase!

# Use np.lookfor() to look for numpy functions
np.lookfor('mean')
np.mean(a)

# Use wildcard with ? to look for function names
np.me*?

# Making arrays
a = np.array(range(4))
a.ndim
a.shape
len(a)

b = np.array([range(3), range(3)])
b.ndim
b.shape
len(b)

np.arange(10)
np.linspace(0, 1, 6)
np.linspace?
# I'm not really sure what "linspace" is supposed to mean, but this gives you a non-int
# range with a 'number of points' argument.
np.linspace(0, 20, 5)
np.arange(20)
# note that the arguments are different than those of the builtin range()

np.ones((3,5,))
np.zeros((2, 3, 4))
np.diag(np.arange(1, 5))
np.eye(3) # array with 1s on diagonal and 0s elsewhere. "eye" like identity matrix (rolling my eyes so bad)

np.random.rand(4)
np.random.randn(4)
# rand() is uniform and randn() is gaussian and who knows why? who will ever remember this?
# I can tell this is a python package written by scientists...

# Be aware of the data type of the array
a = np.arange(10)
a.dtype

a = np.arange(10.0)
a.dtype

# Plotting with matplotlib
x = np.linspace(0, 3, 20)
x.shape
y = np.random.rand(20)
y.shape
plt.plot(x, y)
plt.plot(x, y, 'o')

image = np.random.randn(30, 30)
plt.imshow(image, cmap=plt.cm.hot)
plt.colorbar()
plt.close()

plt.imshow(image, cmap=plt.cm.gray)
plt.colorbar()
plt.close()

# Skipping the stuff on indexing as I mostly remember it from grad school days

## Section 1.3.2
a = np.arange(1,5)
a + 1 # this adds 1 to every element in the array
2 ** a

# array operations are much faster than the pure python equivalent
a = np.arange(10000)
%timeit a + 1 # 16 microsec
l = range(10000)
%timeit [i + 1 for i in l] # 954 microsec

# note that matrix multiplication is a different thing altogether
c = np.ones((3,3))
c * c # this is just array multiplication
c.dot(c) # this is the matrix multiplication

# transcendental functions
a = np.arange(5)
np.sin(a)
np.log(a) # gives a divide-by-0 WARNING where python would give a ZeroDivisionError
5 / 0
np.exp(a)

# shapes need to match sometimes
a = np.arange(4)
a + np.array([1, 2])
# shape (4, ) vs. (2, )

np.triu?
a = np.triu(np.ones((3, 3)), 1) # a 3 x 3 array with 1s in the upper triangle
a.T # transpose; I'm not sure why this doesn't need () for a method call...
np.lookfor('transpose')
# This is a "view". I'm not sure what that means, but it might have something to do with the lack of method call here...

# numpy.linalg has linear algebra functions, BUT scipy.linalg is more likely to be compiled efficiently, so use it there instead

np.allclose?
np.tril?

# Sum by rows or columns
x = np.array(((1, 1,), (2, 2,)))
x.sum(axis=0) # by columns (first dim)
x.sum(axis=1)# by rows (second dim)
x.mean(axis=0)
x.mean(axis=1)
x.min(axis=0)
x.min(axis=1)

np.argmax?
np.lookfor('diff')

# Download an example data file from the lecture notes
import urllib
file_url = 'http://scipy-lectures.org/_downloads/populations.txt'
help(urllib.request.urlretrieve)
urllib.request.urlretrieve(file_url, filename='populations.txt')

# ipython command for viewing the file
!cat populations.txt

data = np.loadtxt('populations.txt')
data
# transpose the data and save each column to its own array
year, hares, lynxes, carrots = data.T

plt.axes([0.2, 0.1, 0.5, 0.8])
plt.plot(year, hares, year, lynxes, year, carrots)
plt.legend(('Hare', 'Lynx', 'Carrot'), loc=(1.05, 0.5))
plt.show()
#close the plot


