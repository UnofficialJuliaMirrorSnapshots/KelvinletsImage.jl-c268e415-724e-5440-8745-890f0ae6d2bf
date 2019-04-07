module KelvinletsImage

    include("triangleInterpolator.jl")
    using Images, ProgressMeter, ImageView, LinearAlgebra
    export KelvinletsObject, grab, scale, pinch, grabRectangle, makeVideo, __applyAreaVariation__, __applyCloserAreaVariation__, __applyPerimeterBasedAreaVariation__, __applyEdgeBasedAreaVariation__, __applyPonderedAvgAreaVariation__
    
    """
        KelvinletsObject(image::AbstractArray{RGB{N0f8}, 2}, ν::Float64, μ::Float64)

    Initializes KelvinletsObject for a given image *image*,
    poisson ratio *ν*
    and elastic shear modulus *μ*

    # Example:
    ```julia-repl
    julia> object = KelvinletsObject(image, 0.4, 1.)
    ```
    """
    struct KelvinletsObject
        a::Float64
        b::Float64
        c::Float64
        sizeX::Int64
        sizeY::Int64
        image::AbstractArray{RGB{N0f8}, 2}
        function KelvinletsObject(image::AbstractArray{RGB{N0f8}, 2},
                                  ν::Float64,
                                  μ::Float64
                )::KelvinletsObject

            a = 1 / (4pi * μ)
            b = a / (4(1 - ν))
            c = 2 / (3a- 2b)

            sizeY, sizeX = size(image)

            new(a, b, c, sizeX, sizeY, image)
        end
    end

    function __applyVariation__(object::KelvinletsObject,
                                variationFunction::Function,
                                retardationFunction::Function,
                                heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)

        for i=1:object.sizeY
            for j=1:object.sizeX
                Δ = variationFunction([i, j])
                
                dx1 = j - 1
                dx2 = object.sizeX - j
                dy1 = i - 1
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy) / object.sizeY
                x = 2(object.sizeX/2 - dx) / object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{N0f8}, 2}
        
        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)

        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)

        for i=1:object.sizeY
            for j=1:object.sizeX

                if i < minY && j < minX #Top-LEFT corner
                    Δ = variationFunction([i, j], [minY, minX])
                elseif i < minY && j >= minX && j <= maxX # Top
                    Δ = variationFunction([i, j], [minY, j])
                elseif i < minY && j > maxX #Top-RIGHT corner
                    Δ = variationFunction([i, j], [minY, maxX])
                elseif j < minX && i >= minY && i <= maxY #LEFT
                    Δ = variationFunction([i, j], [i, minX])
                elseif j >= minX && j <= maxX && i >= minY && i <= maxY # MIDDLE
                    Δ = variationFunction([i, j], [i, j])
                elseif j > maxX && i >= minY && i <= maxY #RIGHT
                    Δ = variationFunction([i, j], [i, maxX])
                elseif i > maxY && j < minX #Lower-LEFT corner
                    Δ = variationFunction([i, j], [maxY, minX])
                elseif i > maxY && j >= minX && j <= maxX #BASE
                    Δ = variationFunction([i, j], [maxY, j])
                elseif i >  maxY && j > maxX #Lower-RIGHT corner
                    Δ = variationFunction([i, j], [maxY, maxX])
                end

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy)/object.sizeY
                x = 2(object.sizeX/2 - dx)/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyCloserAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)

        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)
        for i=1:object.sizeY
            for j=1:object.sizeX

                df00 = norm([i, j] - [minY, minX])
                df01 = norm([i, j] - [minY, maxX])
                df10 = norm([i, j] - [maxY, minX])
                df11 = norm([i, j] - [maxY, maxX])

                closer = min(df00, df01, df10, df11)

                if closer == df00
                    Δ = variationFunction([i, j], [minY, minX])
                elseif closer == df01
                    Δ = variationFunction([i, j], [minY, maxX])
                elseif closer == df10
                    Δ = variationFunction([i, j], [maxY, minX])
                elseif closer == df11
                    Δ = variationFunction([i, j], [maxY, maxX])
                end

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - dy)/object.sizeY
                x = 2(object.sizeX/2 - dx)/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)
                
                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]

                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)

    end

    function __applyPerimeterBasedAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                divided = [0., 0.]
                divisor = 0.

                for k = minX:maxX
                    d = norm([i, j] - [minY, k]) + 1
                    divided += (1/d) * [minY, k]
                    divisor += 1/d
                end

                for k = minX+1:maxX-1
                    d = norm([i, j] - [maxY, k]) + 1
                    divided += (1/d) * [maxY, k]
                    divisor += 1/d
                end

                for k = minY+1:maxY
                    d = norm([i, j] - [k, minX]) + 1
                    divided += (1/d) * [k, minX]
                    divisor += 1/d
                end

                for k = minY+1:maxX
                    d = norm([i, j] - [k, maxX]) + 1
                    divided += (1/d) * [k, maxX]
                    divisor += 1/d
                end

                newReference = divided/divisor
                Δ = variationFunction([i, j], round.(Int64, newReference))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i

                dx = min(dx1, dx2)
                dy = min(dy1, dy2)

                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX

                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)
        
                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyEdgeBasedAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                d1 = norm([i, j] - [minY, minX]) + 1
                d2 = norm([i, j] - [minY, maxX]) + 1
                d3 = norm([i, j] - [maxY, minX]) + 1
                d4 = norm([i, j] - [maxY, maxX]) + 1

                newReference = ((1/d1) * [minY, minX] + (1/d2) * [minY, maxX] + (1/d3) * [maxY, minX] + (1/d4) * [maxY, maxX]) / ((1/d1) + (1/d2) + (1/d3) + (1/d4))
                Δ = variationFunction([i, j], round.(Int64, newReference))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i
        
                dx = min(dx1, dx2)
                dy = min(dy1, dy2)
        
                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX
        
                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)


                maxnorm = norm([object.sizeY, object.sizeX])

                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)


                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __applyPonderedAvgAreaVariation__(object::KelvinletsObject,
                                    points::Array{Int64, 2},
                                    variationFunction::Function,
                                    retardationFunction::Function,
                                    heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        allΔ = zeros(object.sizeY, object.sizeX, 2)
        numOfPixels = (maxY - minY) * (maxX - minX)

        testImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)

        @showprogress for i=1:object.sizeY
            for j=1:object.sizeX

                d1 = norm([i, j] - [minY, minX]) + 1
                d2 = norm([i, j] - [minY, maxX]) + 1
                d3 = norm([i, j] - [maxY, minX]) + 1
                d4 = norm([i, j] - [maxY, maxX]) + 1

                deform1 = variationFunction([i,j], [minY, minX])
                deform2 = variationFunction([i,j], [minY, maxX])
                deform3 = variationFunction([i,j], [maxY, minX])
                deform4 = variationFunction([i,j], [maxY, maxX])

                Δ = ((1/d1) * deform1 + (1/d2) * deform1 + (1/d3) * deform1 + (1/d4) * deform1) /
                     ((1/d1) + (1/d2) + (1/d3) + (1/d4))

                dx1 = j
                dx2 = object.sizeX - j
                dy1 = i
                dy2 = object.sizeY - i
        
                dx = min(dx1, dx2)
                dy = min(dy1, dy2)
        
                y = 2(object.sizeY/2 - (dy - 1))/object.sizeY
                x = 2(object.sizeX/2 - (dx - 1))/object.sizeX
        
                Δ[1] *= retardationFunction(y)
                Δ[2] *= retardationFunction(x)

                maxnorm = norm([object.sizeY, object.sizeX])
                testImg[i, j] = RGB{N0f8}(norm(Δ)/maxnorm, norm(Δ)/maxnorm, norm(Δ)/maxnorm)

                Δ += [i, j]
                allΔ[i, j, 1] = Δ[1]
                allΔ[i, j, 2] = Δ[2]
            end
        end
        if heatmap
            return testImg
        end
        return __interpolateVariation__(object, allΔ)
    end

    function __interpolateVariation__(object::KelvinletsObject,
                                      allΔ::Array{Float64, 3}
            )::Array{RGB{N0f8}, 2}

        interpImg = fill(RGB(0, 0, 0), object.sizeY, object.sizeX)
        
        rasterize = function(A, B, C)
            colorA = object.image[A[1], A[2]]
            colorB = object.image[B[1], B[2]]
            colorC = object.image[C[1], C[2]]
            triangleInterpolator.rasterizationBBOX(interpImg,
                            [allΔ[A[1], A[2], 1], allΔ[A[1], A[2], 2]],
                            [allΔ[B[1], B[2], 1], allΔ[B[1], B[2], 2]],
                            [allΔ[C[1], C[2], 1], allΔ[C[1], C[2], 2]],
                            colorA, colorB, colorC
            )
        end
        
        for i=1:object.sizeY-1
          for j=1:object.sizeX-1
            rasterize([i, j], [i, j+1], [i+1, j+1])
            rasterize([i, j], [i+1, j], [i+1, j+1])
          end
        end
        return interpImg
    end

    """
        grab(object::KelvinletsObject, x0::Array{Int64}, force::Array{Float64}, ϵ::Float64, heatmap::Bool)

    grabs a point in an image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a force vector *force*,
    a brush size *ϵ*
    and whether you would like to generate the heatmap of the deformation
    # Example:
    ```julia-repl
        julia> newImage = grab(obj, [100, 100], [100., 0.], 50.)
    ```
    """    
    function grab(object::KelvinletsObject,
                  x0::Array{Int64},
                  force::Array{Float64},
                  ϵ::Float64,
                  heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        grabFunc = function(x::Array{Int64})
            
            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end
        
        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, grabFunc, retardationFunc, heatmap)
    end

    """
        scale(object::KelvinletsObject, x0::Array{Int64}, scale::Float64, ϵ::Float64, heatmap::Bool)

    Scales a point on an image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a scale alpha *scale*,
    a brush size *ϵ*,
    and whether you would like to generate the heatmap of the deformation
    # Example:
    ```julia-repl
        julia> newImage = scale(obj, [100, 100], 50., 50.)
    ```
    """
    function scale(object::KelvinletsObject,
                   x0::Array{Int64},
                   force::Float64,
                   ϵ::Float64,
                   heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        scaleFunc = function(x::Array{Int64})
            
            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)

            return (2 * object.b - object.a) *
                   ( (1 / rϵ^2) +
                   ((ϵ^2)) / (2 * (rϵ^4))) *
                   (force * r)
        end
        
        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, scaleFunc, retardationFunc, heatmap)
    end
    
    """
         pinch(object::KelvinletsObject, x0::Array{Int64}, force::Array{Float64, 2}, ϵ::Float64, heatmap::Bool)

    Pinches the image given a KelvinletsObject *obj*,
    a pressure point *x0*,
    a force matrix *force*,
    a brush size *ϵ*,
    and whether you would like to generate the heatmap of the deformation
    # Example:
    ```julia-repl
        julia> newImage = pinch(obj, [100, 100], [10000. 0. ; 10000. 0.], 50.)
    ```
    """
    function pinch(object::KelvinletsObject,
                   x0::Array{Int64},
                   force::Array{Float64, 2},
                   ϵ::Float64,
                   heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        pinchFunc = function(x::Array{Int64})
            
            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            return  -2 * object.a * ((1 / rϵ^2) +
                    (ϵ^2 / rϵ^4)) * force * r +
                    4 * object.b * ((1 / rϵ^2) * force -
                    (1 / rϵ^4) * (r' * force * r) * I) * r
        end
        
        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, pinchFunc, retardationFunc, heatmap)
    end

    """
        grabRectangle(object::KelvinletsObject, points::Array{Int64, 2}, force::Array{Float64}, ϵ::Float64, heatmap::Bool)

    Grabs the image using a square brush given a KelvinletsObject *obj*,
    a given matrix representation of a rectangle *points* in the following format --> [minY, minX ; maxY, maxX],
    a force vector *force*,
    a given brush size *ϵ* (to calculate te variation of outside pixels),
    and whether you would like to generate the heatmap of the deformation
    # Example:
    ```julia-repl
        julia> newImage = grabRectangle(obj, [100 100 ; 250 250], [100., 0.], 50.)
    ```
    """
    function grabRectangle(object::KelvinletsObject,
                           points::Array{Int64, 2},
                           force::Array{Float64},
                           ϵ::Float64,
                           areaVariation::Function,
                           heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        grabFunc = function(x::Array{Int64},
                            x0::Array{Int64}
                    )
        
            r = x - x0
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return areaVariation(object, points, grabFunc, retardationFunc, heatmap)

    end

    function __r__(a, b, c, d, x)
        ax = a - x
        bx = b - x
        cx = c - x
        dx = d - x
        
        mainDiag = a - d
               
        Aa = abs(ax[1]) * abs(ax[2])
        Ab = abs(bx[1]) * abs(bx[2])
        Ac = abs(cx[1]) * abs(cx[2])
        Ad = abs(dx[1]) * abs(dx[2])
        
        At = abs(mainDiag[1]) * abs(mainDiag[2])
        
        da = norm(ax) + 1
        db = norm(bx) + 1
        dc = norm(cx) + 1
        dd = norm(dx) + 1
        
        return 1/512 * ((At/(Aa + Ab + Ac + Ad)) - 1) * 
                   (x - 
                    ((1/da) * a + (1/db) * b + (1/dc) * c + (1/dd) * d)) /
                    ((1/da) + (1/db) + (1/dc) + (1/dd)
                   )
    end


    function grabRectangle__new(object::KelvinletsObject,
                           points::Array{Int64, 2},
                           force::Array{Float64},
                           ϵ::Float64,
                           heatmap::Bool
            )::Array{RGB{N0f8}, 2}

        minX, maxX = points[:, 1]
        minY, maxY = points[:, 2]

        a = [minY, minX]
        b = [minY, maxX]
        c = [maxY, minX]
        d = [maxY, maxX]

        grabFunc = function(x::Array{Int64})

            r = __r__(a, b, c, d, x)
            rLength = norm(r)
            rϵ = sqrt(rLength^2 + ϵ^2)
            kelvinState = (((object.a - object.b)/rϵ) * I +
                            (object.b / rϵ^3) * r * r' +
                            (object.a / 2) * (ϵ^2 / rϵ^3) * I)
            object.c * ϵ * kelvinState * force
        end

        retardationFunc = α -> (cos(π * α) + 1) / 2
        return __applyVariation__(object, grabFunc, retardationFunc, heatmap)

    end

    """
        makeVideo(object::KelvinletsObject, kelvinletsFunction::Function, x0::Array{Int64}, force, ϵ::Float64, frames::Int64)

    Computes a video for a given KelvinletsImage deformation function *KelvinletsFunction*, 
    using a *KelvinletsObject* as reference,
    a *x0* point of force application,
    a *force* matrix/vector (depending on the function),
    a *ϵ* brush size
    and a number of frames *frames*
    # Example
    ```julia-repl
        julia> video = (object, grab, [100, 100], [200., 0.], 70., 20)
    ```
    """
    function makeVideo(object::KelvinletsObject,
                       kelvinletsFunction::Function,
                       x0::Array{Int64},
                       force,
                       ϵ::Float64,
                       frames::Int64
        )::Array{RGB{N0f8},3}
        
        if typeof(force) == Float64
            var = range(0, stop=force, length=frames)
        else    
            var = range(fill(0, size(force)), stop=force, length=frames)    
        end
        
        video = Array{RGB{N0f8}}(undef, object.sizeY, object.sizeX, frames)
        @showprogress for i=1:frames
            video[:,:,i] = kelvinletsFunction(object, x0, var[i], ϵ, false)
        end
        imshow(video)
        return video
    end
end
