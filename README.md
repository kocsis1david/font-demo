# Lorem ipsum font demo

This repository contains a Vulkan demo of the font renderer I made, an explanation of how it works can be found [here](howitworks.md).

![Screenshot](images/screenshot.png)

## Building the project

### Windows

Download and set the correct include and library directories to the dependencies: FreeType, GLFW, Vulkan SDK.

Run the demo passing the path to a font as the first argument e.g. "C:\windows\fonts\times.ttf"

### Linux

Assuming FreeType, GLFW and the Vulkan SDK are installed.

Use CMake and Make to build the project:

```sh
cmake -S . -B build
make -C build
```

Run the demo passing in a font of your choice:

```sh
cd build
./font-demo "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"
```

### macOS

This has been tested on macOS 10.14.2 with version 1.1.92.1 of the LunarG Vulkan SDK.

Install CMake and Freetype:

```sh
brew install cmake freetype
```

Install glfw **from HEAD** because version 3.2 of glfw does not include Vulkan support for macOS:

```sh
brew install glfw --HEAD
```

Download the latest [LunarG Vulkan SDK](https://vulkan.lunarg.com/sdk/home) and set up the environment variables as per the instructions [in the docs](https://vulkan.lunarg.com/doc/view/latest/mac/getting_started.html#user-content-command-line).

Use CMake and Make to build the project:

```sh
cmake -S . -B build
make -C build
```

Run the demo passing in a font of your choice:

```sh
cd build
./font-demo "/Library/Fonts/Times New Roman.ttf" 
```
