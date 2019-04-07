---
Author: Guilherme Gomes Haetinger
Title: KelvinletsImage.jl
Date: Nov. 8th, 2018
---
# KelvinletsImage.jl

Implementation for Kelvinlets Deformations presented on 
[Regularized Kelvinlets: Sculpting Brushes based on Fundamental Solutions of Elasticity](https://graphics.pixar.com/library/Kelvinlets/paper.pdf) 
from *Fernando De Goes* and *Doug L. James* on **Julia v1.0**.

Deformations for the *Grab*, *Scale* and *Pinch* brushes. Including a **Rectangle** brush for the *Grab* function.

---

## Usage

### Initialization

You must first initialize the structure **KelvinletsObject** to start deforming a given image. To do that you must first select an image end define a *poisson ratio* **ν** and an *elastic shear modulus* **μ** (both of them are explained on the above paper).


```julia
    using KelvinletsImage, TestImages

    ν = 0.4
    μ = 1.0
    image = testimages("mandrill")

    object = KelvinletsObject(image, ν, μ)
```

![original image](./sampleImages/original.png)

After initializing the presented object, you will use ir for further operations.

### Grab Brush


```julia
    using KelvinletsImage

    pressurePosition = [256, 256]
    forceVector = [200., 0.]
    ϵ = 70. # Brush Size
    grabbedImage = grab(object, pressurePosition, forceVector, ϵ)
```

![grabbed image](./sampleImages/grab.png)

```julia
    using KelvinletsImage

    frames = 20
    grabbedImageGIF = makeVideo(object, grab, pressurePosition, forceVector, ϵ, frames)
```

![grabbed image GIF](./sampleImages/grabVid.gif)


### Scale Brush


```julia
    using KelvinletsImage

    pressurePosition = [256, 256]
    scale = -200000. # Negative value = inflates .. Positive Value = Contracts
    ϵ = 70. # Brush Size
    scaledImage = scale(object, pressurePosition, forceVector, ϵ)
```

![grabbed image](./sampleImages/scale.png)

```julia
    using KelvinletsImage

    frames = 20
    grabbedImageGIF = makeVideo(object, scale, pressurePosition, scale, ϵ, frames)
```

![grabbed image GIF](./sampleImages/scaleVid.gif)


### Pinch Brush


```julia
    using KelvinletsImage

    pressurePosition = [256, 256]
    forceMatrix = [0. 0.; 0. 300000.]
    ϵ = 300. # Brush Size
    grabbedImage = pinch(object, pressurePosition, forceMatrix, ϵ)
```

![grabbed image](./sampleImages/pinch.png)

```julia
    using KelvinletsImage

    frames = 20
    grabbedImageGIF = makeVideo(object, pinch, pressurePosition, forceMatrix, ϵ, frames)
```

![grabbed image GIF](./sampleImages/pinchVid.gif)


### Grab Rectangle Brush


```julia
    using KelvinletsImage

    pressureRectangle = [156 156; 256 256]
    forceVector = [150., 150.]
    ϵ = 50. # Brush Size
    grabbedImage = pinch(object, pressureRectangle, forceVector, ϵ)
```

![grabbed image](./sampleImages/grabRectangle.png)

```julia
    using KelvinletsImage

    frames = 20
    grabbedImageGIF = makeVideo(object, grabRectangle, pressureRectangle, forceVector, ϵ, frames)
```

![grabbed image GIF](./sampleImages/grabRectangle.gif)
