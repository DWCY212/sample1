import rpy2.robjects as ro

# Run R code
ro.r('''
x <- rnorm(100)
summary(x)
''')

# Capture R output and print
print(ro.r['summary'](ro.r['x']))
