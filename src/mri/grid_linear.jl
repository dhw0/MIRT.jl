#=
    Project:      MIRTdev
    File:         grid_linear.jl
    Description:  Crude "gridding" based on linear interpolation
    Author:       Daniel Wan
    Date Created: 2019-10-01
    Date Updated: 2019-10-21
=#

export mri_grid_linear

using Interpolations
using FFTW
using NFFT
using ImageView
using Colors
using Gtk.ShortNames
using LinearAlgebra
include("mri_objects.jl")

"""
`function  mri_grid_linear(kspace, ydata, N, fov)
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
# Copyright 2019-10-01; Jeff Fessler; University of Michigan

function mri_grid_linear(kspace, ydata, N, fov)

    length(N) == 1 && (N = [N N]) 
    length(N) == 2 || error("bad N")
    length(fov) == 1 && (fov = [fov fov])
    length(fov) == 2 || error("bad fov")

    kg = [(-N[1] / 2 : N[1] / 2 - 1) / fov[1]]
    push!(kg, (-N[2] / 2 : N[2] / 2 - 1) / fov[2])
    # need to put into one array http://juliamath.github.io/Interpolations.jl/latest/control/#Gridded-interpolation-1
    knots = (kspace[:,1], kspace[:,2])
    itp = interpolate(knots, ydata, Gridded(Linear()))
    etp = extrapolate(itp, 0)
    yhat = [etp(i, j) for i in kg[1], j in kg[2]]
  
    xg = [(-N[1] / 2 : N[1] / 2 - 1) / N[1] * fov[1]]
    push!(xg, (-N[2] / 2 : N[2] / 2 - 1) / N[2] * fov[2])
    dk = 1 ./ fov
    xhat = fftshift(ifft(fftshift(yhat))) * prod(dk) * prod(N)
    tmp1 = sinc.(xg[1] / fov[1])
    tmp2 = sinc.(xg[2] / fov[2])
    xhat = xhat ./ (tmp1 * tmp2'); # gridding post-correction for linear interp
    return xhat, yhat, xg, kg
end
    
    #
    # test function
    #
function mri_grid_linear(test::Symbol)
    if (test != :test)
        error("argument must be \":test\" if only one argument is specified")
    end

    fov = 256 # [mm] typical brain FOV
    N0 = 64 # nominal image size()
    
    t = range(0, stop = N0 / 2 * 2 * pi, length = N0 ^ 2) # crude spiral:
    kspace = N0 / 2 * (1 / fov) * [cos.(t) sin.(t)] .* (t[:, 1] / maximum(t))
    
    Ndisp = 256; # display images like this...
    x1d = (-Ndisp / 2 : Ndisp / 2 - 1) / Ndisp * fov
    x1dd = [i for i in x1d, j in x1d]
    x2dd = [j for i in x1d, j in x1d]

    # TODO: find what mri_objects does
    #=
    obj = mri_objects["case1"]
    xtrue = obj.image(x1dd, x2dd)
    ytrue = obj.kspace[kspace[:,1], kspace[:,2]]
    
    im plc 2 3
    clim = [0 2]
    im[1, x1d, x1d, xtrue, "x true", clim], cbar
    
    Ng = 128
    (xhat, yhat, xg, kg) = mri_grid_linear(kspace, ytrue, Ng, fov)
    
    im[2, kg[1], kg[2], abs(yhat), "|y_{grid]|"], cbar
    im[3, xg[1], xg[2], abs(xhat), "|x| \"gridding\"", clim], cbar
    =#

    img1 = Gray.(rmul!(xtrue, 1 / maximum(xtrue)))
    img2 = Gray.(rmul!(abs(yhat), 1 / maximum(abs(yhat))))
    img3 = Gray.(rmul!(abs(xhat), 1 / maximum(abs(xhat))))

    grid, frames, canvases = canvasgrid((1,3))  # 1 row, 3 columns
    ImageView.imshow(canvases[1,1], img1)
    imshow(canvases[1,2], img2)
    imshow(canvases[1,3], img3)
    
    win = Window(grid)
    Gtk.showall(win)
end