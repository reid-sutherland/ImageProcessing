# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /private/var/folders/zr/ftj_dq_50t1795b82w5hgvhm0000gn/T/AppTranslocation/58EF6B35-40B1-416A-83CB-1C68F44168D6/d/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /private/var/folders/zr/ftj_dq_50t1795b82w5hgvhm0000gn/T/AppTranslocation/58EF6B35-40B1-416A-83CB-1C68F44168D6/d/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Reid/ImageProcessing/Assignment2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Reid/ImageProcessing/Assignment2/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Assignment2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Assignment2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Assignment2.dir/flags.make

CMakeFiles/Assignment2.dir/main.cpp.o: CMakeFiles/Assignment2.dir/flags.make
CMakeFiles/Assignment2.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Reid/ImageProcessing/Assignment2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Assignment2.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Assignment2.dir/main.cpp.o -c /Users/Reid/ImageProcessing/Assignment2/main.cpp

CMakeFiles/Assignment2.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Assignment2.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Reid/ImageProcessing/Assignment2/main.cpp > CMakeFiles/Assignment2.dir/main.cpp.i

CMakeFiles/Assignment2.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Assignment2.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Reid/ImageProcessing/Assignment2/main.cpp -o CMakeFiles/Assignment2.dir/main.cpp.s

# Object files for target Assignment2
Assignment2_OBJECTS = \
"CMakeFiles/Assignment2.dir/main.cpp.o"

# External object files for target Assignment2
Assignment2_EXTERNAL_OBJECTS =

Assignment2: CMakeFiles/Assignment2.dir/main.cpp.o
Assignment2: CMakeFiles/Assignment2.dir/build.make
Assignment2: /usr/local/lib/libopencv_stitching.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_superres.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_videostab.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_aruco.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_bgsegm.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_bioinspired.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_ccalib.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_dnn_objdetect.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_dpm.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_face.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_fuzzy.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_hfs.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_img_hash.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_line_descriptor.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_optflow.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_reg.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_rgbd.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_saliency.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_stereo.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_structured_light.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_surface_matching.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_tracking.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_xfeatures2d.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_ximgproc.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_xobjdetect.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_xphoto.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_shape.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_photo.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_dnn.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_datasets.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_ml.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_plot.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_video.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_calib3d.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_features2d.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_highgui.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_videoio.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_phase_unwrapping.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_flann.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_imgcodecs.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_objdetect.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_imgproc.3.4.1.dylib
Assignment2: /usr/local/lib/libopencv_core.3.4.1.dylib
Assignment2: CMakeFiles/Assignment2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Reid/ImageProcessing/Assignment2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Assignment2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Assignment2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Assignment2.dir/build: Assignment2

.PHONY : CMakeFiles/Assignment2.dir/build

CMakeFiles/Assignment2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Assignment2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Assignment2.dir/clean

CMakeFiles/Assignment2.dir/depend:
	cd /Users/Reid/ImageProcessing/Assignment2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Reid/ImageProcessing/Assignment2 /Users/Reid/ImageProcessing/Assignment2 /Users/Reid/ImageProcessing/Assignment2/cmake-build-debug /Users/Reid/ImageProcessing/Assignment2/cmake-build-debug /Users/Reid/ImageProcessing/Assignment2/cmake-build-debug/CMakeFiles/Assignment2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Assignment2.dir/depend

