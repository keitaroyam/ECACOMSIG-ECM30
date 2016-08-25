"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License.
"""

import numpy
from cctbx import miller
from cctbx import crystal

class Predictions:
    def __init__(self):
        pass

    def read_xparm(self, xpin):
        fin = open(xpin)
        assert "XPARM.XDS" in fin.readline() # Skip first header line

        # Line 2
        sp = fin.readline().split() # starting frame & angle, osc. range, rotation axis
        self.starting_frame, self.starting_angle, self.osc_range = int(sp[0]), float(sp[1]), float(sp[2])
        m2 = numpy.array(map(float, sp[3:6])) # rotation axis

        # Line 3
        sp = fin.readline().split() # wavelength & incident beam direction
        self.wavelength = float(sp[0])
        incident_beam = numpy.array(map(float, sp[1:4]))
        self.s0 = incident_beam / numpy.linalg.norm(incident_beam) / self.wavelength # |S0| = 1/lambda

        m = numpy.empty(dtype=numpy.float, shape=(3,3)) # (m1 m2 m3) matrix
        m[:,1] = m2 / numpy.linalg.norm(m2)
        m[:,0] = numpy.cross(m[:,1], self.s0) # m1 = m2 x s0 / |m2 x s0|
        m[:,0] /= numpy.linalg.norm(m[:,0])
        m[:,2] = numpy.cross(m[:,0], m[:,1]) # m3 = m1 x m2
        self.m_matrix = m
        
        # Line 4
        sp = fin.readline().split() # Space group number & unit cell constants
        self.crystal_symmetry = crystal.symmetry(unit_cell=map(float, sp[1:7]),
                                                 space_group_symbol=int(sp[0]))

        # Line 5,6,7 real space vectors
        a_axis = map(float, fin.readline().split())
        b_axis = map(float, fin.readline().split())
        c_axis = map(float, fin.readline().split())
        
        b_mat = numpy.array([a_axis, b_axis, c_axis]).transpose() # (b0 b1 b2) matrix
        self.astar_matrix = numpy.linalg.inv(b_mat) # (b0* b1* b2*)^t matrix

        # Line 8 detector dimensions and pixel size
        sp = fin.readline().split()
        self.nxy = map(int, sp[1:3]) # NX, NY
        self.qxy = map(float, sp[3:5]) # QX, QY

        # Line 9 ORGX,Y & detector distance
        sp = fin.readline().split()
        self.orgxy = map(float, sp[:2])
        self.detector_F = float(sp[2])
        
        # Line 10,11,12
        d1 = map(float, fin.readline().split())
        d2 = map(float, fin.readline().split())
        d3 = map(float, fin.readline().split()) # this is actually d1 x d2

        self.d_matrix = numpy.array([d1,d2,d3]).transpose()  # (d1 d2 d3) matrix
    # read_xparm()

    def prep_indices(self, d_min, d_max=None):
        mset = miller.build_set(self.crystal_symmetry, anomalous_flag=True, d_min=d_min, d_max=d_max)
        return mset.indices()
    # prep_indices()

    def calc_centroids(self, indices):
        h = numpy.array(indices) # hkl in each row
        m, d, a, F = self.m_matrix, self.d_matrix, self.astar_matrix, self.detector_F
        x0, y0 = self.orgxy
        qx, qy = self.qxy
        s0 = self.s0

        p0s = numpy.dot(h, a) # p0* vector in each row
        p0s_m = numpy.dot(p0s, m) # p0* with m1,m2,m3 basis

        s0_m = numpy.dot(m.transpose(), s0) # s0 with m1,m2,m3 basis

        p0s_lensq = numpy.sum(p0s**2, axis=1) # |p0*|^2 = |p*|^2

        ps_m = numpy.empty(p0s_m.shape) # p* with m1,m2,m3 basis
        ps_m[:,2] = (-0.5*p0s_lensq - s0_m[1]*p0s_m[:,1]) / s0_m[2]
        ps_m[:,1] = p0s_m[:,1]
        ps_m[:,0] = p0s_lensq - ps_m[:,1]**2 - ps_m[:,2]**2 # sqrt after check

        sel_ok = ps_m[:,0] > 0 # No solution (blind region) if < 0
        h, p0s_m, ps_m = h[sel_ok], p0s_m[sel_ok], ps_m[sel_ok]
        ps_m[:,0] = numpy.sqrt(ps_m[:,0])
        
        self.predicted_hkl = numpy.empty((0, 3), dtype=numpy.int) # h,k,l
        self.predicted_data = numpy.empty((0, 4)) # x, y, phi, zeta

        for p1sign in (+1, -1):
            ps_m[:,0] *= p1sign
            phi = numpy.arctan2(p0s_m[:,2]*ps_m[:,0] - p0s_m[:,0]*ps_m[:,2], # sin(phi) rho^2
                                p0s_m[:,0]*ps_m[:,0] + p0s_m[:,2]*ps_m[:,2]) # cos(phi) rho^2
            e1_m = numpy.cross(ps_m, s0_m)
            e1_m /= numpy.linalg.norm(e1_m, axis=1).reshape(e1_m.shape[0], 1)
            zeta = e1_m[:,1] # as e1_m is expressed in m1,m2,m3 system (zeta = m2 . e1)

            s = s0 + numpy.dot(ps_m, m.transpose())
            s_d = numpy.dot(s, d) # S vector with d1,d2,d3 basis
            sel_ok = F*s_d[:,2] > 0
            s_d, h_ok, phi, zeta = s_d[sel_ok], h[sel_ok], phi[sel_ok], zeta[sel_ok]
            xdet = x0 + F*s_d[:,0]/s_d[:,2]/qx
            ydet = y0 + F*s_d[:,1]/s_d[:,2]/qy

            self.predicted_hkl = numpy.row_stack([self.predicted_hkl, h_ok])
            self.predicted_data = numpy.row_stack([self.predicted_data, numpy.column_stack([xdet, ydet, phi, zeta])])

    # calc_centroids()

    def get_predicted_positions(self, sigma_m, frame, esd_factor=3):
        phi = self.starting_angle + self.osc_range * (frame - self.starting_frame + 0.5) # Is this correct?
        print "  Phi at frame %d = %.3f" % (frame, phi)
        
        phi, sigma_m, osc_range = numpy.deg2rad([phi, sigma_m, self.osc_range])

        phi_calc = self.predicted_data[:,2]
        zeta = self.predicted_data[:,3]

        phi_diff = numpy.fmod(phi_calc - phi, 2.*numpy.pi)
        phi_diff[phi_diff < -numpy.pi] += 2.*numpy.pi
        phi_diff[phi_diff > numpy.pi] -= 2.*numpy.pi

        sel = numpy.abs(phi_diff) < osc_range/2. + esd_factor * sigma_m / numpy.abs(zeta)

        return self.predicted_hkl[sel], self.predicted_data[sel]
    # get_predicted_positions()

# class Predictions

if __name__ == "__main__":
    import sys

    xparm_in = sys.argv[1]
    d_min = float(sys.argv[2])
    sigma_m = float(sys.argv[3])
    frames = map(int, sys.argv[4:])

    
    preds = Predictions()
    print "Reading XPARM.XDS.."
    preds.read_xparm(xparm_in)

    indices = preds.prep_indices(d_min)

    print "Calculating the predicted centroids.."
    preds.calc_centroids(indices)
    
    for frame in frames:
        print "Calculating the predictions on frame %d.." % frame
        pindices, pdata = preds.get_predicted_positions(sigma_m, frame)

        ofs = open("prediction_%.6d.adx" % frame, "w")
        for (h,k,l), (x, y, phi, zeta) in zip(pindices, pdata):
            ofs.write("%d %d %d %d %d\n"%(x,y, h,k,l))

        ofs.close()
        
