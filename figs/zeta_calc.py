import numpy

d1 = numpy.array([1.,0.,0.])
d2 = numpy.array([0.,1.,0.])
d3 = numpy.array([0.,0.,1.])
s0 = numpy.array([0.,0.,1.])
m2 = numpy.array([1.,0.,0.])

print "x y zeta"
for x in xrange(-100, 101):
  for y in xrange(-100, 101):
    s = x*d1 + y*d2 + 100*d3
    e1 = numpy.cross(s, s0)
    e1 /= numpy.linalg.norm(e1)
    zeta = abs(numpy.dot(e1, m2))
    print x,y,zeta
