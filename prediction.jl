#=
This is just directly translated from prediction.py.
Not optimized for Julia (actually a bit slower than NumPy code).
- As advised by Takanori, by not using mapslices and lambda expressions, calc_centroids() became about x2 faster (now a bit faster than NumPy code)
=#

using PyCall

type Xparm
    starting_frame::Int
    starting_angle::Float64
    osc_range::Float64
    wavelength::Float64
    s0::Array{Float64,1}
    m_matrix::Array{Float64,2}
    unit_cell::Array{Float64,1}
    sgnum::Int
    astar_matrix::Array{Float64,2}
    qxy::Array{Float64,1}
    nxy::Array{Float64,1}
    orgxy::Array{Float64,1}
    detector_F::Float64
    d_matrix::Array{Float64,2}
end

function prep_indices(unit_cell, sgnum, d_min, d_max=nothing)
    @pyimport cctbx.miller as miller
    @pyimport cctbx.crystal as crystal
    xs = crystal.symmetry(tuple(unit_cell...), sgnum)
    mset = miller.build_set(xs, anomalous_flag=true, d_min=d_min, d_max=d_max)
    mi = mset["indices"]()
    ret = zeros(Int, mi["size"](), 3)
    for i in 1:size(ret, 1)
        for j in 1:3
            ret[i, j] = mi[i][j]
        end
    end
    return ret
end

function read_xparm(xpin)
    open(xpin) do fin
        @assert contains(readline(fin), "XPARM.XDS") # Skip first header line

        # Line 2
        sp = split(readline(fin)) # starting frame & angle, osc. range, rotation axis
        starting_frame, starting_angle, osc_range = parse(Int, sp[1]), float(sp[2]), float(sp[3])
        m2 = map(float, sp[4:6]) # rotation axis
        
        # Line 3
        sp = split(readline(fin)) # wavelength & incident beam direction
        wavelength = float(sp[1])
        incident_beam = map(float, sp[2:4])
        s0 = incident_beam / norm(incident_beam) / wavelength # |S0| = 1/lambda
        
        m = zeros(Float64, 3, 3) # (m1 m2 m3) matrix
        m[1:3, 2] = m2 / norm(m2)
        m[1:3, 1] = cross(m[1:3, 2], s0) # m1 = m2 x s0 / |m2 x s0|
        m[1:3, 1] /= norm(m[1:3, 1])
        m[1:3, 3] = cross(m[1:3, 1], m[1:3, 2]) # m3 = m1 x m2
        m_matrix = m
        
        # Line 4
        sp = split(readline(fin)) # Space group number & unit cell constants
        unit_cell = map(float, sp[2:7])
        sgnum = parse(Int, sp[1])

        # Line 5,6,7 real space vectors
        a_axis = map(float, split(readline(fin)))
        b_axis = map(float, split(readline(fin)))
        c_axis = map(float, split(readline(fin)))
        
        b_mat = [a_axis b_axis c_axis] # (b0 b1 b2) matrix
        astar_matrix = inv(b_mat) # (b0* b1* b2*)^t matrix

        # Line 8 detector dimensions and pixel size
        sp = split(readline(fin))
        nxy = map(x->parse(Int,x), sp[2:3]) # NX, NY
        qxy = map(float, sp[4:5]) # QX, QY

        # Line 9 ORGX,Y & detector distance
        sp = split(readline(fin))
        orgxy = map(float, sp[1:2])
        detector_F = float(sp[3])
        
        # Line 10,11,12
        d1 = map(float, split(readline(fin)))
        d2 = map(float, split(readline(fin)))
        d3 = map(float, split(readline(fin))) # this is actually d1 x d2

        d_matrix = [d1 d2 d3]  # (d1 d2 d3) matrix

        return Xparm(starting_frame, starting_angle, osc_range,
                     wavelength, s0, m_matrix, unit_cell, sgnum,
                     astar_matrix, qxy, nxy, orgxy, detector_F, d_matrix)
    end
end

function calc_centroids(xparm, h)
    # h = (hkl) in each row
    m, d, a, F = xparm.m_matrix, xparm.d_matrix, xparm.astar_matrix, xparm.detector_F
    x0, y0 = xparm.orgxy
    qx, qy = xparm.qxy
    s0 = xparm.s0

    p0s = h * a # p0* vector in each row
    p0s_m = p0s * m # p0* with m1,m2,m3 basis

    s0_m = transpose(m) * s0 # s0 with m1,m2,m3 basis

    p0s_lensq = sum(p0s.^2, 2) # |p0*|^2 = |p*|^2 #

    ps_m = zeros(size(p0s_m)) # p* with m1,m2,m3 basis
    ps_m[:,3] = (-0.5*p0s_lensq - s0_m[2]*p0s_m[:,2]) / s0_m[3]
    ps_m[:,2] = p0s_m[:,2]
    ps_m[:,1] = p0s_lensq - ps_m[:,2].^2 - ps_m[:,3].^2 # sqrt after check
    
    sel_ok = ps_m[:,1] .> 0 # No solution (blind region) if < 0 # HOW!?
    h, p0s_m, ps_m = h[sel_ok,:], p0s_m[sel_ok,:], ps_m[sel_ok,:]
    ps_m[:,1] = sqrt(ps_m[:,1])
    
    predicted_hkl = Array(Int, (0,3)) # h,k,l
    predicted_data = Array(Float64, (0,4)) # x, y, phi, zeta

    for p1sign in (+1, -1)
        ps_m[:,1] *= p1sign
        phi = atan2(p0s_m[:,3].*ps_m[:,1] - p0s_m[:,1].*ps_m[:,3], # sin(phi) rho^2
                    p0s_m[:,1].*ps_m[:,1] + p0s_m[:,3].*ps_m[:,3]) # cos(phi) rho^2

        e1_m = zeros(Float64, size(ps_m))
        for i in 1:size(ps_m, 1)
            e1_m[i,:] = cross(ps_m[i,:][:], s0_m)
        end
        
            
        e1_m = e1_m ./ sqrt(sum(e1_m.^2, 2))
        zeta = e1_m[:,2] # as e1_m is expressed in m1,m2,m3 system (zeta = m2 . e1)

        s = ps_m * transpose(m)
        for i in 1:size(ps_m, 1)
            s[i,:] += s0'
        end
        #
        s_d = s * d # S vector with d1,d2,d3 basis
        sel_ok = F*s_d[:,3] .> 0
        s_d, h_ok, phi, zeta = s_d[sel_ok,:], h[sel_ok,:], phi[sel_ok,:], zeta[sel_ok,:]
        xdet = x0 + F*s_d[:,1] ./ s_d[:,3]/qx
        ydet = y0 + F*s_d[:,2] ./ s_d[:,3]/qy

        predicted_hkl = vcat(predicted_hkl, h_ok)
        predicted_data = vcat(predicted_data, [xdet ydet phi zeta])
    end
    return predicted_hkl, predicted_data
end

function get_predicted_positions(xparm, predicted_hkl, predicted_data, sigma_m, frame, esd_factor=3)
    phi = xparm.starting_angle + xparm.osc_range * (frame - xparm.starting_frame + 0.5) # Is this correct?
    println(@sprintf("  Phi at frame %d = %.3f", frame, phi))
    
    phi, sigma_m, osc_range = deg2rad([phi, sigma_m, xparm.osc_range])

    phi_calc = predicted_data[:,3]
    zeta = predicted_data[:,4]

    phi_diff = mod(phi_calc - phi, 2pi)
    phi_diff[phi_diff .< -pi] += 2pi
    phi_diff[phi_diff .> pi] -= 2pi

    sel = abs(phi_diff) .< osc_range/2. + esd_factor * sigma_m ./ abs(zeta)

    return predicted_hkl[sel,:], predicted_data[sel,:]
end # get_predicted_positions()


function main()
    xparm_in = ARGS[1]
    d_min = float(ARGS[2])
    sigma_m = float(ARGS[3])
    frames = map(x->parse(Int, x), ARGS[4:length(ARGS)])

    xp = read_xparm(xparm_in)

    indices = prep_indices(xp.unit_cell, xp.sgnum, d_min)

    println("Calculating the predicted centroids..")
    #@time calc_centroids(xp, indices)
    #@time calc_centroids(xp, indices)
    predicted_hkl, predicted_data = calc_centroids(xp, indices)
    
    for frame in frames
        println(@sprintf("Calculating the predictions on frame %d..", frame))
        pindices, pdata = get_predicted_positions(xp, predicted_hkl, predicted_data, sigma_m, frame)

        open(@sprintf("prediction_%06d.adx", frame), "w") do ofs
            for i in 1:size(pindices, 1)
                hkl = pindices[i,:]
                data = pdata[i,:]
                write(ofs, @sprintf("%.1f %.1f %d %d %d\n", data[1:2]..., hkl...))
            end
        end
    end
end

main()
