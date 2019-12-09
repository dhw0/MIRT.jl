using SpecialFunctions #necessary to grab Bessel functions
"""
function jinc(x)
returns jinc(x), defined as J1(x)/x, where J1 is a Bessel function of the first kind.
No particular options.
"""
function jinc(x::Real) #calculates jinc(x)
    if(x == 0) return pi/4 end
    y = abs(x)
    return besselj1(pi*y) / (2*y)
end
