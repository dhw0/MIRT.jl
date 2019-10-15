export mri_objects
#note: how do we use fftw
using FFTW
#we haven't defined st in here, so we'll need that obj def for the code to compile.

#changes here: we define an array of states, different structs that hold data. 

#structs here
#this might not be an efficient way to do this?
struct circ2
    data
end
struct dirac2
    data
end
struct gauss2
    data
end
struct rect2
    data
end
struct cycl3
    data
end
struct dirac3
    data
end
struct gauss3
    data
end
struct rect3
    data
end

#NOTE: NEED TO IMPLEMENT ERROR AND OTHER STUFF
function mri_objects(x...)
    if size(x) < 1
        #replaces help(filename), might need more refinement
        throw(DomainError(x,"input must have options argument"));
    end
    if size(x) == 1 && x[1] == "test"
        mri_objects_test()
        return;
    end
    while(size(x))
        #looks like a recursive call thru the program.
        #note that we'll be filling the data array here
    end
    if(size(x) % 2)
        throw(DomainError(x,"bad args. remainder supplied."));
    end    
end
#end todo

#the main functions we're dealing with
function mri_objects_kspace(data,u,v,w)
    out = zeros(size(u))
    for i = 1:size(data)
        out += mri_objects_kspace_typerun(data[i],u,v,w)
    end
    return out
end
function mri_objects_image(data,x,y,z) #if we need 2d, just make the third coord
    out = zeros(size(x))
    for i = 1:size(data) #data is 1d because it's an array of data structs
        out += mri_objects_kspace_typerun(data[i],x,y,z)
    end
    return out
end
#back compatability in case ppl want to run 2d scans. 
#These will work on partially 3d things as well; they'll just drop the additional one.
function mri_objects_kspace(data,u,v)
    return mri_objects_kspace(data,u,v,0)
end
function mri_objects_image(data,x,y)
    return mri_objects_image(data,x,y,0)
end


#The big chunk of functions for each struct type.
function mri_objects_image_typerun(data::gauss2,x,y,z)
    if(size(data.data,2) != 5)
        throw(DomainError(data,"gauss2 requires 5 parameters"));
    end
    z = zeros(size(data.data,1))
    zw = z + Base.Inf
    data.data = [data.data[:,1:2] z data.data[:,3:4] zw data.data[:,5]]
    convert(gauss3,data) #convert to an arg that will run dirac3
    return mri_objects_image_typerun(data,x,y,z)
end
function mri_objects_image_typerun(data::gauss3,x,y,z)
    out = zeros(size(x));
    if(size(data.data,2) != 7)
        throw(DomainError(data,"gauss3 requires 7 parameters"));
    end
    for i=1:size(data.data,1)
        par = data.data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4] / sqrt(log(256)) * sqrt(2*pi)
        yw = par[5] / sqrt(log(256)) * sqrt(2*pi)
        zw = par[6] / sqrt(log(256)) * sqrt(2*pi)
        out = out + par[7] .* exp(-pi * ((x-xc)/xw).^2) .* exp(-pi * ((y-yc)/yw).^2) .* exp(-pi * ((z-zc)/zw).^2);
    end
    return out
end
function mri_objects_image_typerun(data::circ2,x,y,z)
    if(size(data.data,2) != 4)
        throw(DomainError(data,"circ2 requires 4 parameters"));
    end
    data.data = [data.data[:,1:2] z data.data[:,3] 1+z data.data[:,4]]
    convert(cycl3,data) #convert to an arg that will run dirac3
    return mri_objects_image_typerun(data,x,y,z)
end
function mri_objects_image_typerun(data::cycl3,x,y,z)
    out = zeros(size(x))
    if(size(data.data,2) != 6)
        throw(DomainError(data,"cycl3 requires 6 parameters"));
    end
    for i = 1:size(data.data,1)
        par = data.data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        rad = par[4]
        len = par[5]
        #what is rect? Need to replace.
        #none of this stuff is in the matlab documentation lmao
        out = out + par(6) * ((x-xc).^2+(y-yc).^2 < rad^2) .* rect((z-zc)/len);
    end
    return out
end
function mri_objects_image_typerun(data::rect2,x,y,z)
    if(size(data.data,2) != 5)
        throw(DomainError(data,"rect2 requires 5 parameters"));
    end
    data.data = [data.data[:,1:2] z data.data[:,3:4] 1+z data.data[:,5]]
    convert(rect3,data) #convert to an arg that will run dirac3
    return mri_objects_image_typerun(data,x,y,z)
end
function mri_objects_image_typerun(data::rect3,x,y,z)
    for i=1:size(data.data,1)
        par = data.data[i,:]
        xc = par[1]
        yc = par[2]
        zc = par[3]
        xw = par[4]
        yw = par[5]
        zw = par[6]
        #again, we need rect to exist
        out = out + par[7] .* rect((x-xc)/xw) .* rect((y-yc)/yw) .* rect((z-zc)/zw);
    end
end
function mri_objects_image_typerun(data::dirac2,x,y,z)
    if(size(data.data,2) != 3)
        throw(DomainError(data,"dirac2 requires 3 parameters"));
    end
    data.data = [data.data[:,1:2] zeros(size(data.data,0),1) data.data[:,3]]
    convert(dirac3,data) #convert to an arg that will run dirac3
    return mri_objects_image_typerun(data,x,y,z)
end

function mri_objects_image_typerun(data::dirac3,x,y,z)
    if(size(data.data,2) != 4)
        throw(DomainError(data,"dirac3 requires 4 parameters"));
    end
    out = zeros(size(x))
    for i=1:size(data.data,1)
        out += data.data[i,4]*((x == data.data[i,1]) .* (y == data.data[i,2]) .* (z == data.data[i,3]))
    end
    out(out != 0) = Base.Inf
    return out
end
#TODO: what does out(out != 0) = inf, warn do? Ask fess.
function mri_objects_kspace_typerun(data::gauss2,u,v,w)
    if(size(data.data,2) != 5)
        throw(DomainError(data,"gauss2 requires 5 parameters"));
    end
    z = zeros(size(data.data,1));
    zw = (1+z) * sqrt(log(256)) / sqrt(2*pi);
    data.data = [data.data[:,1:2] z data.data[:,3:4] zw data.data[:,5]]
    convert(gauss3,data) #convert to an arg that will run gauss3
    return mri_objects_kspace_typerun(data,u,v,w)
end
function mri_objects_kspace_typerun(data::gauss3,u,v,w)
    
end
function mri_objects_kspace_typerun(data::circ2,u,v,w)
    if(size(data.data,2) != 4)
        throw(DomainError(data,"circ2 requires 4 parameters"));
    end
    z = zeros(size(data:data,1))
    data.data = [data.data[:,1:2] z data.data[:,3] 1+z data.data[:,4]]
    convert(cycl3,data) #convert to an arg that will run cycl3
    return mri_objects_kspace_typerun(data,u,v,w)
end
function mri_objects_kspace_typerun(data::cycl3,u,v,w)
    out = zeros(size(u))
    if(size(data.data,2) != 6)
        throw(DomainError(data,"cycl3 requires 6 parameters"));
    end
    for i = 1:size(data.data,1)
        par = data.data[i,:]
        xc = par[1];
        yc = par[2];
        zc = par[3];
        rad = par[4];
        len = par[5];
        #what are jinc and nufft_sinc?
        out = out + par(6) * rad^2 * 4*jinc(2*sqrt(u.^2+v.^2)*rad) .* (len * nufft_sinc(w*len)) .* exp(-2i*pi*(u*xc + v*yc + w*zc));
    end
    return out
end
function mri_objects_kspace_typerun(data::rect2,u,v,w)
    
end
function mri_objects_kspace_typerun(data::rect3,u,v,w)
    
end
function mri_objects_kspace_typerun(data::dirac2,u,v,w)
    if(size(data.data,2) != 3)
        throw(DomainError(data,"dirac2 requires 3 parameters"));
    end
    data.data = [data.data[:,1:2] zeros(size(data.data,0),1) data.data[:,3]]
    convert(dirac3,data) #convert to an arg that will run dirac3
    return mri_objects_kspace_typerun(data,u,v,w)
end
function mri_objects_kspace_typerun(data::dirac3,u,v,w)
    if(size(data.data,2) != 4)
        throw(DomainError(data,"dirac3 requires 4 parameters"));
    end
    out = zeros(size(u))
    for i=1:size(data.data,1)
        out += data[i,4] * Base.MathConstants.e ^ (-2im*pi*(u*data.data[i,1] + v*data.data[i,2] + w*data.data[i,3]))
    end
    return out
end
function mri_objects_case1(fov,arg)
    #rect2
    rp = [
	[0 0		200 200	1];
	[-50 -50		40 40	1];
	[50 -50		20 20	1];
	[0 50		50 50	1];
];
#gauss2
gp = [
	[-70 0		1 1 1];
	[-60 0		2 2 1];
	[-50 0		3 3 1];
	[-40 0		4 4 1];
	[-20 0		5 5 1];
	[00 0		6 6 1];
	[20 0		7 7 1];
	[50 0		8 8 1];
];
#NOTE: No easy access to isVar. Bring up with Fessler.
if(arg == "cm")
    rp(:,1:6) /= 10
    gp(:,1:6) /= 10
end
return mri_objects('rect2',rp,'gauss2',gp)
end

function mri_objects_test4(fov,arg)
    #cycl3
 cp = [0 0 0 fov(1)*.4 fov(3)*1 1] #not sure why the fov's multiplied by 1. Bring up.
 #rect3
 rp = [
    [-50 -50   0	40 40 40	1];
    [50 -50  40	20 20 50	1];
    [0  50 -40	30 30 60	1];
 ]
 #I'm willing to replace these, but make sure that this works. Bring up.
 rp(:,1:3) = rp(:,1:3)/256 .* [fov fov fov]
 rp(:,4:6) = rp(:,4:6)/256 .* [fov fov fov]

gp = [
	[-70 0 0		1 1 1	1];
	[-60 0 0		2 2 2	1];
	[-50 0 0		3 3 3	1];
	[-40 0 0		4 4 4	1];
	[-20 0 0		5 5 5	1];
	[00 0 0		6 6 6	1];
	[20 0 0		7 7 7	1];
	[50 0 0		8 8 8	1];
];
#again, bring these up.
gp(:,1:3) = gp(:,1:3)/256 .* [fov fov fov fov fov fov fov fov] #ha
gp(:,4:6) = gp(:,4:6)/256 .* [fov fov fov fov fov fov fov fov]
#again, need to check what isvar is doing here and the best way to replicate it. 
if(arg == "cm")
    rp(:,1:6) /= 10
    gp(:,1:6) /= 10
end
return ['cycl3',cp,'rect3',rp,'gauss3',gp]
end

function mri_objects_test2()
#note: need to implement image_geom
ig = image_geom('nx',2^7,'offsets','dst','dx',4)
if 1 #why??????? ASK Fess
    shift = [0.1 0.2] .* ig.fovs;
    sizes = [.15 .1] .* ig.fovs;
    tests = ['circ2',[shift sizes[1] 2];'gauss2',[shift sizes 2];'rect2',[shift sizes 2]]
    #WHAT IS A PLC
    for i = 1:size(tests,1)
        otype = tests[i,1]
        param = tests[i,2]
        st = mri_objects(otype,param)
        i3 = mri_objects_image(st,ig.xg,ig.yg,ig.zg)
        s3 = 
        s3 *= abs(ig.dx * ig.dy * ig.dz)
        fg = ig.fg
        f2 = mri_objects_kspace(st,fg[:])
        #MISSING FUNCTION MAX_PERCENT_DIFF

        #ALSO TODO: graphing

    end
end
#something printing-related here. TODO: find prompt replacements, among other things.
st = mri_objects('case1')
xt = mri_objects_image(st,ig.xg,ig.yg)
#again, what's im. Need to get graphing down.
end
#2d tests
function mri_objects_test3()
    #again, need to define image_geom
    ig = image_geom('nx', 2^7, 'nz', 2^5, 'offsets', 'dsp', 'dx', 4, 'dz', 3);
    if 1 #why
        shift = [0.1 0.2 0.3] .* ig.fovs;
        sizes = [0.15 0.1 0.2] .* ig.fovs;
        tests = ['cyl3', [shift sizes([1 3]) 2];
            'gauss3', [shift sizes 2];
            'rect3', [shift sizes 2];
        ]
        #again, img
        for i = 1 : size(tests,1)
            otype = tests[i,1];
            param = tests[i,2];
            st = mri_objects(otype,param)
            i3 = mri_objects_image(st,ig.xg,ig.yg,ig.zg)
            s3 = 
            s3 *= abs(ig.dx * ig.dy * ig.dz)
            fg = ig.fg
            f3 = mri_objects_kspace(st,fg[:])
            #missing function mpd, as before
            
            #also todo: graphing

        end
    end

    st = mri_objects('fov',ig.fovs,'test4')
    xt = mri_objects_image(st,ig.xg,ig.yg,ig.zg)
    
end

function mri_objects_test
    mri_objects_test2()
    mri_objects_test3()
end
