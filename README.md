# GeneralBrokenLines

![GBL logo](GblLogo.png)

A trajectory based on General Broken Lines is a track refit to add the description of multiple scattering to an initial trajectory based on the propagation in a magnetic field (and average energy loss). It is constructed from a sequence of (pairs of) thin scatterers describing the multiple scattering in the material between adjacent measurement planes.

Detailed documentation can be found at https://millepede.pages.desy.de/general-broken-lines/


## Obtaining the code 

We recommend using the latest tagged release - currently version **04-00-03**. 
This can be done via 
```
wget https://gitlab.desy.de/millepede/general-broken-lines/-/archive/V04-00-03/general-broken-lines-V04-00-03.tar.gz
```

Alternatively, you can clone the sources from the package's [gitlab](https://gitlab.desy.de/millepede/general-broken-lines) repository or the [GitHub mirror](https://github.com/GeneralBrokenLines/GeneralBrokenLines). 

GBL is also distributed with in the [SPack](https://spack.io/) package manager. 

## Building the C++ package

### Requirements

To build the code, you need to have CMake > 3.14 and Eigen3.
Optionally, ROOT is supported. 

The package also relies on the [Mille](https://gitlab.desy.de/millepede/mille) library to communicate with Millepede-II - this is automatically downloaded and installed in case no existing installation is found. 

### Compiling the package 

The most frequently used version of GBL is the C++ implementation in the `cpp` directory. 
This can be built using the `CMake` toolkit: 

1) Create a compilation directory, e.g. `build` and change into it
    `mkdir build; cd build`
2) Configure the build using the `CMake` command. Additional flags can be used to specify the installation location and fine-tune the build. See the `CMake` documentation for more details. An option specific to GLB is that by passing `-DSUPPORT_ROOT=on`, you can enable support for ROOT in GBL. 
    `cmake ../general-broken-lines/cpp`
3) start the build process
    `make && make install`

The installation folder will, if not manually specified using the `CMAKE_INSTALL_PREFIX` flag, default to a subfolder `GBLInstall` within the build directory.

### Activating the environment 

To make your system pick up the library and the example programs, you can call the `gblsetup.sh` script generated inside your install directory. 

### Building the documentation 

To build the documentation, you need doxygen (version > 1.7.0 recommended ) on your system.
Invoke, also in the build directory:
  `make doc`

## Contributing

Contributions in form of issues or merge requests are very welcome! 
If you do not have access to the DESY GitHub instance, please feel free to use the [GitHub mirror](https://github.com/GeneralBrokenLines/GeneralBrokenLines) to submit a PR. 
