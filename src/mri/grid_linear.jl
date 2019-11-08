#=
    Project:      MIRTdev
    File:         grid_linear.jl
    Description:  Crude "gridding" based on linear interpolation
    Author:       Daniel Wan
    Date Created: 2019-10-01
    Date Updated: 2019-10-21
=#
# Copyright 2019-10-01; Jeff Fessler; University of Michigan


export mri_grid_linear

using Interpolations
using FFTW
using NFFT
using LinearAlgebra
using MIRT:jim
using Plots
using Dierckx
include("mri_objects.jl")

"""
    (xhat, yhat, xg, kg) = mri_grid_linear(kspace, ydata, N, fov)

very crude "gridding" based on linear interpolation.
not recommended as anything but a straw man | perhaps
for initializing iterative methods.
in
    kspace		[M 2]	kx & ky k-space sample locations
    ydata		[M]	complex Fourier sample values
    N			desired image size()
    fov			field of view()
		kspace & fov must be in reciprocal units!
out
    xhat		[N N]	image estimate
    yhat		[N N]	gridding k-space data
    xg	([N1], [N2])	object-space uniform sample locations
    kg	([N1], [N2])	k-space uniform sample locations in 1d
"""
function mri_grid_linear(kspace, ydata, N, fov)

    length(N) == 1 && (N = [N N]) 
    length(N) == 2 || throw("bad N")
    length(fov) == 1 && (fov = [fov fov])
    length(fov) == 2 || throw("bad fov")

    kg = [(-N[1] / 2 : N[1] / 2 - 1) / fov[1], 
          (-N[2] / 2 : N[2] / 2 - 1) / fov[2]]
    # need to put into one array http://juliamath.github.io/Interpolations.jl/latest/control/#Gridded-interpolation-1
    # figure out how to interpolate scattered data
    
    yreal = real.(ydata)
    yimag = imag.(ydata)
    yreal = [i for i in yreal, j in yreal]
    knots = (kspace[:,1], kspace[:,2])
    sort(knots[1])
    sort(knots[2])

    itpreal = interpolate(knots, yreal, Gridded(Linear()))
    itpreal = extrapolate(itpreal, 0)
    yhatr = [itpreal(i, j) for i in kg[1], j in kg[2]]

    

    spl = Spline2D(kspace[:,1], kspace[:,2], ydata, kx = 1, ky = 1, s = 0.0)
    yhat = evalgrid(spl, kg[1], kg[2])
    xg = [(-N[1] / 2 : N[1] / 2 - 1) / N[1] * fov[1], 
          (-N[2] / 2 : N[2] / 2 - 1) / N[2] * fov[2]]
    dk = 1 ./ fov
    xhat = fftshift(ifft(fftshift(yhat))) * prod(dk) * prod(N)
    tmp1 = sinc.(xg[1] / fov[1])
    tmp2 = sinc.(xg[2] / fov[2])
    xhat = xhat ./ (tmp1 * tmp2') # gridding post-correction for linear interp
    return xhat, yhat, xg, kg
end
    

"""
    function  mri_grid_linear(test::Symbol)

Test function for mri_grid_linear.
Displays three images of uniform data from non-uniform data.
in
    test        Symbol  symbol to indicate test routine
out
    N/A
"""
function mri_grid_linear(test::Symbol)
    if (test != :test)
        error("argument must be \":test\" if only one argument is specified")
    end

    fov = 256 # [mm] typical brain FOV
    N0 = 64 # nominal image size()
    
    t = range(0, stop = N0 / 2 * 2 * pi, length = N0 ^ 2) # crude spiral:
    kspace = N0 / 2 * (1 / fov) * [cos.(t) sin.(t)] .* (t[:, 1] / maximum(t))
    
    Ndisp = 256 # display images like this...
    x1d = (-Ndisp / 2 : Ndisp / 2 - 1) / Ndisp * fov
    x1dd = [i for i in x1d, j in x1d]
    x2dd = [j for i in x1d, j in x1d]

    #obj = mri_objects("case1")
    xtrue = mri_objects_image(mri_objects(:fov, [0], :case1), x1dd, x2dd)
    ytrue = mri_objects_kspace(mri_objects(:fov, [0], :case1), kspace[:,1], kspace[:,2])
    
    Ng = 128
    (xhat, yhat, xg, kg) = mri_grid_linear(kspace, ytrue, Ng, fov)
    
    p1 = jim(xtrue, x = x1d, y = x1d, clim = (0, 2), title = "x true")
    p2 = jim(abs(yhat), x = kg[1], y = kg[2], title = "|y_{grid]|")
    p3 = jim(abs(xhat), x = xg[1], y = xg[2], clim = (0, 2), title = "|x| \"gridding\"")
    p = plot(p1, p2, p3)
end