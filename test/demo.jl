
# Pionier: K = 2190nm ± 390nm
# Gravity: H = 1630nm ± 307nm
# Matisse: L = 3450nm ± 472nm
#          M = 4750nm ± 460nm
#          N = 10500nm ± 2500nm

using EasyFITS, FITSHeaders, Unitful
using Revise, ImageMetrics
using ImageMetrics.InterferometricImagingContest2024: ContestData, resample_cube

CONTEST_DIR = "/home/eric/data/interferometry_imaging_beauty_contest/2024/"
data = Dict{String,Any}()
for key in (
    "ContestData/Obj1/output_Spiral_1710269539282_H",
    #"ContestData/Obj1/output_Spiral_1710273340532_K",
    #"ContestData/Obj1/output_Spiral_1710273706182_L",
    #"ContestData/Obj1/output_Spiral_1710274957235_N",
    #"ContestData/Obj2/Model_H",
    #"ContestData/Obj2/Model_K",
    #"ContestData/Obj2/Model_L",
    #"ContestData/Obj2/Model_N",
    "Drevon/OBJ1/POLY_OBJ1_GRAV_cube",
    #"Drevon/OBJ1/POLY_OBJ1_MATLM_cube",
    #"Drevon/OBJ1/POLY_OBJ1_MATN_cube",
    #"Drevon/OBJ1/POLY_OBJ1_PION_cube",
    #"Drevon/OBJ2/POLY_OBJ2_GRAV_cube",
    #"Drevon/OBJ2/POLY_OBJ2_MATLM_cube",
    #"Drevon/OBJ2/POLY_OBJ2_MATN_cube",
    #"Drevon/OBJ2/POLY_OBJ2_PION_cube",
    "Norris/obj1_gravity_norris",
    #"Norris/obj1_matisse29_norris",
    #"Norris/obj1_matisse8_norris",
    #"Norris/obj1_pionier_norris",
    #"Norris/obj2_gravity_norris",
    #"Norris/obj2_matisse3_norris",
    #"Norris/obj2_matisse8_norris",
    #"Norris/obj2_pionier_norris",)
    data[key] = ContestData(joinpath(CONTEST_DIR, key*".fits"))
end


A = data["ContestData/Obj1/output_Spiral_1710269539282_H"];
B = data["Drevon/OBJ1/POLY_OBJ1_GRAV_cube"];
C = data["Norris/obj1_gravity_norris"];

Bmax = 70u"m"
resolution(λ::Unitful.Length) = uconvert(u"rad", λ/3Bmax)

B_img = resample_cube(B.data;
                      input_pixelsize = B.pixelsize,
                      output_pixelsize = A.pixelsize,
                      input_wave = B.wave,
                      output_wave = B.wave,
                      output_resolution = resolution.(B.wave));
B_ref = resample_cube(A.data;
                      input_pixelsize = A.pixelsize,
                      output_pixelsize = A.pixelsize,
                      input_wave = A.wave,
                      output_wave = B.wave,
                      output_resolution = resolution.(B.wave))

function load_norris_data(filename::AbstractString, )
    openfits(joinpath(NORRIS_DIR, filename)) do file
        hdu = file[1]
        naxis1 = hdu["NAXIS1"].integer
        naxis2 = hdu["NAXIS2"].integer
        cdelt1 = hdu["CDELT1"].float
        cdelt2 = hdu["CDELT2"].float
        @assert abs(cdelt1) ≈ abs(cdelt2)
        A = ImageMetrics.InterferometricImagingContest2024.resample_cube(
        read(Array{Float64}, hdu);
        input_pixelsize = abs(cdelt1),
        input_wave = ((1:naxis3) .- crpix3).*cdelt3 .+ crval3,
        output_pixelsize = ,
                       output_wave::AbstractVector{<:Number},
        output_resolution
h = FitsHeader(hdu)
